#ifdef GINT_HAVE_LIBINT
#include "libint.hh"

namespace gint {
namespace gaussian {
namespace libint {

namespace detail {
std::vector<libint2::Atom> make_libint_atomlist(const Molecule& molecule) {
  std::vector<libint2::Atom> atom_list;
  atom_list.reserve(molecule.n_atoms());

  for (const gint::Atom& atom : molecule) {
    libint2::Atom a_converted;
    a_converted.x = atom.x;
    a_converted.y = atom.y;
    a_converted.z = atom.z;
    a_converted.atomic_number = static_cast<int>(atom.nuclear_charge);

    atom_list.push_back(std::move(a_converted));
  }

  return atom_list;
}

std::vector<std::pair<double, std::array<double, 3>>> inline make_libint_point_charges(
      const Molecule& molecule) {
  std::vector<std::pair<double, std::array<double, 3>>> ret;
  ret.reserve(molecule.n_atoms());
  for (const auto& atom : molecule) {
    ret.emplace_back(static_cast<double>(atom.nuclear_charge),
                     std::array<double, 3>{{atom.x, atom.y, atom.z}});
  }
  return ret;
}
}  // namespace detail

//
// LibintSystem
//

LibintSystem::LibintSystem(const std::string& basisset_name,
                           krims::SubscriptionPointer<const Molecule> structure_ptr)
      : m_basis(basisset_name, detail::make_libint_atomlist(*structure_ptr)),
        m_structure_ptr(std::move(structure_ptr)),
        m_point_charges(detail::make_libint_point_charges(*m_structure_ptr)) {}

//
// LibintShellData
//
krims::Range<size_t> LibintShellData::shell_range(
      krims::Range<size_t> bfct_indices) const {
  // We assume here that shell2bf is sorted:
  assert_dbg(std::is_sorted(std::begin(shell2bf), std::end(shell2bf)),
             krims::ExcInternalError());
  assert_dbg(!shell2bf.empty(), krims::ExcInternalError());

  // Get the lower bound of the element one greater than the basis function index.
  // The lower bound is the first element which is not smaller than index+1,
  // so the one before that is the first shell we care about.
  const auto lower = std::lower_bound(std::begin(shell2bf), std::end(shell2bf),
                                      bfct_indices.lower_bound() + 1);
  const size_t begin_shell =
        lower == std::end(shell2bf)
              ? shell2bf.size() - 1
              : static_cast<size_t>(lower - std::begin(shell2bf)) - 1;

  // Get the lower bound on the upper bound of the basis function indices.
  // This is the first shell we do not care about any more
  const auto upper =
        std::lower_bound(lower, std::end(shell2bf), bfct_indices.upper_bound());
  const size_t end_shell = upper == std::end(shell2bf)
                                 ? shell2bf.size()
                                 : static_cast<size_t>(upper - std::begin(shell2bf));

  return {begin_shell, end_shell};
}

//
// LibintCollection
//

const std::string IntegralCollection::id = "gaussian/libint";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : m_system{}, m_global{} {
  // Try to setup the basis of the system to model.
  const std::string& basis_set = parameters.at<const std::string>("basis_set");
  auto structure_ptr = parameters.at_ptr<const Molecule>("structure");

  const std::string error_msg =
        "Could not construct gaussian basis for the given structure and the "
        "gaussian basis set name \"" +
        basis_set + "\". Check that all atoms are supported for a basis of this name.";
  const std::string details = "\n\nDetails:\n--------\n";

  try {
    m_system = LibintSystem(basis_set, structure_ptr);
  } catch (std::exception& e) {
    assert_throw(false, ExcInvalidIntegralParameters(error_msg + details + e.what()));
  } catch (std::string& e) {
    assert_throw(false, ExcInvalidIntegralParameters(error_msg + details + e));
  } catch (...) {
    assert_throw(false, ExcInvalidIntegralParameters(error_msg));
  }
}

Integral<real_stored_mtx_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  using libint2::Operator;

  switch (type) {
    case IntegralType::nuclear_attraction:
      return make_integral<OneElecIntegralCore>(Operator::nuclear, m_system, m_global);
    case IntegralType::overlap:
      return make_integral<OneElecIntegralCore>(Operator::overlap, m_system, m_global);
    case IntegralType::kinetic:
      return make_integral<OneElecIntegralCore>(Operator::kinetic, m_system, m_global);
    case IntegralType::coulomb:
      return make_integral<ERICore>(/* exchange= */ false, m_system, m_global);
    case IntegralType::exchange:
      return make_integral<ERICore>(/* exchange= */ true, m_system, m_global);
  }

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<real_stored_mtx_type>(nullptr);
}

//
// OneElecIntegralCore
//

template <typename Kernel>
void OneElecIntegralCore::compute_kernel(const krims::Range<size_t>& rows,
                                         const krims::Range<size_t>& cols,
                                         Kernel&& kernel) const {
  // TODO libint2::Engine is not thread-safe, so for running this in
  // parallel, we need one engine object for each thread.
  //
  // TODO Parallelise here!

  const libint2::BasisSet& basis = base_type::system().basis();
  const size_t max_nprim = basis.max_nprim();         // Maximum number of primitives
  const int max_l = static_cast<int>(basis.max_l());  // Maximum angular momentum
  const int derivative_order = 0;                     // Calculate no derivatives
  const real_type tolerance = std::numeric_limits<real_type>::epsilon();
  libint2::Engine int_engine(m_operator, max_nprim, max_l, derivative_order, tolerance);

  // If this is a nuclear attraction operator, set the point charges:
  if (m_operator == libint2::Operator::nuclear) {
    int_engine.set_params(base_type::system().point_charges());
  }

  // Loop over shell sets {s,t}
  for (size_t s : rows) {
    for (size_t t : cols) {
      // Compute the integral for a shell pair and call the kernel function with the
      // location where the computed integrals are stored.
      int_engine.compute(basis[s], basis[t]);
      kernel(s, t, int_engine.results()[0]);
    }  // t
  }    // s
}

scalar_type OneElecIntegralCore::operator()(size_t row, size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  const LibintShellData data{base_type::system()};
  scalar_type ret = 0;
  auto kernel = [&data, &row, &col, &ret](size_t s, size_t t, const real_type* values) {
    assert_dbg(data.first_bfct(s) <= row, krims::ExcInternalError());
    assert_dbg(data.first_bfct(t) <= col, krims::ExcInternalError());

    if (values == nullptr) return;  // All values are zero => ret will be zero, too

    // Find the one value we need and extract it.
    const size_t fs = row - data.first_bfct(s);
    const size_t ft = col - data.first_bfct(t);
    ret = values[fs * data.n_bfct(t) + ft];
  };

  // Compute the one shell pair of the kernel we truly need:
  compute_kernel(data.shell_range({row, row + 1}), data.shell_range({col, col + 1}),
                 std::move(kernel));
  return ret;
}

void OneElecIntegralCore::extract_block(stored_matrix_type& M, const size_t start_row,
                                        const size_t start_col,
                                        const linalgwrap::Transposed mode,
                                        const scalar_type c_A,
                                        const scalar_type c_M) const {
  assert_finite(c_A);
  assert_finite(c_M);
  assert_greater_equal(start_row + M.n_rows(), n_rows());
  assert_greater_equal(start_col + M.n_cols(), n_cols());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  // Set elements of M to zero (if c_M == 0)
  // or scale them according to c_M.
  // This deals entirely with the coefficient c_M
  if (c_M == 0) {
    M.set_zero();
  } else {
    M *= c_M;
  }
  if (c_A == 0) return;  // If c_A == 0 we are done

  // The range of basis function indices we are interested in:
  const auto row_range = krims::range(start_row, start_row + M.n_rows());
  const auto col_range = krims::range(start_col, start_col + M.n_rows());

  const LibintShellData data{base_type::system()};
  auto kernel = [&data, &row_range, &col_range, &M, &c_A](size_t s, size_t t,
                                                          const real_type* values) {
    if (values == nullptr) return;  // The whole shell-pair was screened away

    // The first basis function indices we care about in this shell pair block:
    size_t i = std::max(data.first_bfct(s), row_range.lower_bound());
    size_t j = std::max(data.first_bfct(t), col_range.lower_bound());

    // The last basis function indices we a care about in this shell pair block
    const size_t i_end =
          std::min(data.first_bfct(s) + data.n_bfct(s), row_range.upper_bound());
    const size_t j_end =
          std::min(data.first_bfct(t) + data.n_bfct(t), col_range.upper_bound());

    // libint2 packs the integrals into the shellsets in row-major
    // (C-like) form, so this iterates over the results and copy to the matrix M
    for (; i < i_end; ++i) {
      for (; j < j_end; ++j) {
        // Compute offset into values array:
        const size_t fs = i - data.first_bfct(s);
        const size_t ft = j - data.first_bfct(t);
        const size_t fst = fs * data.n_bfct(t) + ft;

        // Shift i and j for placement into M:
        const size_t i_shifted = i - row_range.lower_bound();
        const size_t j_shifted = j - col_range.lower_bound();
        M(i_shifted, j_shifted) = values[fst];
      }  // j
    }    // i

    for (size_t fs = 0; fs < data.n_bfct(s); ++fs) {
      for (size_t ft = 0; ft < data.n_bfct(t); ++ft) {
        const size_t fst = fs * data.n_bfct(t) + ft;  // Index of function pair in result
        const size_t i = data.first_bfct(s) + fs;     // Absolute index in first shell
        const size_t j = data.first_bfct(t) + ft;     // Absolute index in second shell
        M(i, j) = c_A * values[fst];
      }  // j
    }    // i
  };
  compute_kernel(data.shell_range(row_range), data.shell_range(col_range),
                 std::move(kernel));
}

void OneElecIntegralCore::apply(const const_multivector_type& x, multivector_type& y,
                                const linalgwrap::Transposed mode, const scalar_type c_A,
                                const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same since we are symmetric and real, so no
  // switching over mode.

  // Scale the current values of out or set them to zero
  // (if c_y == 0): We are now done with c_y and do not
  // need to worry about it any more in this function
  if (c_y == 0) {
    y.set_zero();
  } else {
    y *= c_y;
  }
  if (c_A == 0) return;  // If c_A == 0 we are done

  const LibintShellData data{base_type::system()};
  auto kernel = [&data, &x, &y, &c_A](size_t s, size_t t, const real_type* values) {
    if (values == nullptr) return;  // The whole shell-pair was screened away

    // libint2 packs the integrals into the shellsets in row-major
    // (C-like) form, so this iterates over the results one after the other
    // and applies it to the vectors
    for (size_t fs = 0; fs < data.n_bfct(s); ++fs) {
      for (size_t ft = 0; ft < data.n_bfct(t); ++ft) {
        const size_t fst = fs * data.n_bfct(t) + ft;  // Index of function pair in result
        const size_t i = data.first_bfct(s) + fs;     // Absolute index in first shell
        const size_t j = data.first_bfct(t) + ft;     // Absolute index in second shell

        for (size_t vec = 0; vec < x.n_cols(); ++vec) {
          y(i, vec) += c_A * values[fst] * x(j, vec);
        }  // vec
      }    // j
    }      // i
  };
  compute_kernel(std::move(kernel));
}

//
// ERI Integral core
//
scalar_type ERICore::operator()(size_t a, size_t b) const {
  assert_greater(a, n_rows());
  assert_greater(b, n_cols());
  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  // TODO   to get it to compile:
  return 0;
}

void ERICore::extract_block(stored_matrix_type& M, const size_t start_row,
                            const size_t start_col, const linalgwrap::Transposed mode,
                            const scalar_type c_A, const scalar_type c_M) const {
  assert_finite(c_A);
  assert_finite(c_M);
  assert_greater_equal(start_row + M.n_rows(), n_rows());
  assert_greater_equal(start_col + M.n_cols(), n_cols());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  // Set elements of M to zero (if c_M == 0)
  // or scale them according to c_M.
  // This deals entirely with the coefficient c_M
  if (c_M == 0) {
    M.set_zero();
  } else {
    M *= c_M;
  }
  if (c_A == 0) return;  // If c_A == 0 we are done

  // TODO
}

void ERICore::apply(const const_multivector_type& x, multivector_type& y,
                    const linalgwrap::Transposed mode, const scalar_type c_A,
                    const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  // Scale the current values of out or set them to zero
  // (if c_y == 0): We are now done with c_y and do not
  // need to worry about it any more in this function
  if (c_y == 0) {
    y.set_zero();
  } else {
    y *= c_y;
  }
  if (c_A == 0) return;  // If c_A == 0 we are done

  // TODO
}

void ERICore::update(const krims::GenMap& map) {
  const std::string occ_coeff_key = Integral<stored_mtx_type>::update_key_coefficients;
  if (!map.exists(occ_coeff_key)) return;

  // Get coefficients as a shared pointer (having ownership)
  coefficients_occupied_ptr =
        static_cast<coefficients_ptr_type>(map.at_ptr<coefficients_type>(occ_coeff_key));

  // We will contract the coefficient row index over the number of
  // basis functions.
  if (coefficients_occupied_ptr->n_vectors() == 0) return;
  assert_size(coefficients_occupied_ptr->n_elem(), base_type::n_bas());
}

}  // namespace libint
}  // namespace gaussian
}  // namespace gint
#endif  // GINT_HAVE_LIBINT
