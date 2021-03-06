//
// Copyright (C) 2017 by the gint authors
//
// This file is part of gint.
//
// gint is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gint is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gint. If not, see <http://www.gnu.org/licenses/>.
//

#include "libint.hh"

#ifdef GINT_HAVE_LIBINT
#include "BasisSet.hh"
#include "IntegralLookupKeys.hh"
#include "Shell.hh"
#include "gint/IntegralUpdateKeys.hh"
#include "read_basisset.hh"
#include <iterator>
#include <krims/DataFiles/FindDataFile.hh>
#include <krims/FileSystem.hh>

namespace gint {
namespace gaussian {
namespace libint {

// TODO This code has some similarities with the code for linking to libcint.
//      Perhaps some more stuff could be generalised and some repeated code could be
//      removed with some work.

namespace detail {
std::vector<std::pair<double, std::array<double, 3>>> inline make_libint_point_charges(
      const Structure& structure) {
  std::vector<std::pair<double, std::array<double, 3>>> ret;
  ret.reserve(structure.n_atoms());
  for (const auto& atom : structure) ret.emplace_back(atom.nuclear_charge, atom.coords);
  return ret;
}

// If boost has the container::small_vector available libint uses this datastructure.
//  Unfortunately a std::vector cannot be directly converted to this datastructure
//  and we need to distinguish between two codepaths, namely one with explict conversion
//  and one without
/** Convert a gint::gaussian::Basis to a libint basis if the vector types agree */
template <bool vectors_types_agree =
                std::is_same<libint2::svector<real_type>, std::vector<real_type>>::value>
std::vector<libint2::Shell> make_libint_basis(Basis basis) {
  std::vector<libint2::Shell> ret;
  ret.reserve(basis.size());

  for (const Shell& sh : basis) {
    libint2::Shell::Contraction cntr{sh.l, sh.pure, std::move(sh.coefficients)};
    libint2::Shell tosh{std::move(sh.exponents), {std::move(cntr)}, std::move(sh.origin)};
    ret.push_back(std::move(tosh));
  }
  return ret;
}

/** Convert a gint::gaussian::Basis to a libint basis if the vector types do not agree */
template <>
std::vector<libint2::Shell> make_libint_basis<false>(Basis basis) {
  std::vector<libint2::Shell> ret;
  ret.reserve(basis.size());

  for (const Shell& sh : basis) {
    libint2::svector<real_type> l_coefficients{
          std::make_move_iterator(sh.coefficients.begin()),
          std::make_move_iterator(sh.coefficients.end())};
    libint2::svector<real_type> l_exponents{std::make_move_iterator(sh.exponents.begin()),
                                            std::make_move_iterator(sh.exponents.end())};

    libint2::Shell::Contraction cntr{sh.l, sh.pure, std::move(l_coefficients)};
    libint2::Shell tosh{std::move(l_exponents), {std::move(cntr)}, std::move(sh.origin)};
    ret.push_back(std::move(tosh));
  }
  return ret;
}

}  // namespace detail

//
// LibintBasisShellData
//

/** When initialised with a LibintSystem it computes and makes available some data
 *  which is needed to work with libint shell pairs
 *
 *  TODO This exists in a similar way in libcint. Perhaps one could create a
 *       generalisation which applies to both.
 */
struct LibintBasisShellData {
  /** Return the index of the first basis function of a particular shell */
  size_t first_bfct(size_t shell) const { return shell2bf[shell]; }

  /** Return the number of basis functions of a particular shell */
  size_t n_bfct(size_t shell) const { return m_system_ptr->basis()[shell].size(); }

  /** Return the range of basis function indices covered by a particular shell */
  krims::Range<size_t> range_bfct(size_t shell) const {
    return {first_bfct(shell), first_bfct(shell) + n_bfct(shell)};
  }

  /** Return the shell information for a particular shell */
  LibintShell shell_info(size_t shell) const {
    return LibintShell{shell, m_system_ptr->basis()[shell].size(), shell2bf[shell]};
  }

  /** Return the number of shells */
  size_t n_shells() const { return m_system_ptr->n_shells(); }

  /** Return the shell range which covers the requested index range
   * in the basis functions.
   *
   * In other words if only some basis function indices are required,
   * this function returns the range of shells which need to be computed
   * in order to compute the values for all basis functions.
   * Note that it is very likely that extra values at the beginning and
   * the end will be computed as well.
   */
  krims::Range<size_t> containing_shell_range(krims::Range<size_t> bfct_indices) const;

  /** The map shell index to index of first basis function of the shell
   *  Obtain the map shell index to basis function index
   *  E.g. shell2bf[0] = index of the first basis function in shell 0
   *       shell2bf[1] = index of the first basis function in shell 1
   */
  const std::vector<size_t> shell2bf;

  explicit LibintBasisShellData(const LibintSystem& system)
        : shell2bf(libint2::BasisSet::compute_shell2bf(system.basis())),
          m_system_ptr("LibintBasisData", system) {}

 private:
  // Pointer to the molecular system info we use:
  krims::SubscriptionPointer<const LibintSystem> m_system_ptr;
};

krims::Range<size_t> LibintBasisShellData::containing_shell_range(
      krims::Range<size_t> bfct_indices) const {
  // We assume here that shell2bf is sorted:
  assert_internal(std::is_sorted(std::begin(shell2bf), std::end(shell2bf)));
  assert_internal(!shell2bf.empty());

  // Shortcut for the full range:
  if (bfct_indices.lower_bound() == 0 &&
      bfct_indices.upper_bound() == static_cast<size_t>(m_system_ptr->n_bas())) {
    return {0, n_shells()};
  }

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
// LibintSystem
//

LibintSystem::LibintSystem(krims::SubscriptionPointer<const Structure> structure_ptr,
                           Basis basis)
      : m_n_bas(0),
        m_basis(detail::make_libint_basis(std::move(basis))),
        m_max_nprim(libint2::BasisSet::max_nprim(m_basis)),
        m_max_l(libint2::BasisSet::max_l(m_basis)),
        m_structure_ptr(std::move(structure_ptr)),
        m_point_charges(detail::make_libint_point_charges(*m_structure_ptr)) {
  m_n_bas = libint2::BasisSet::nbf(m_basis);
}

//
// ERITensor
//

void ERITensor::compute_kernel(const std::array<krims::Range<size_t>, 4>& block,
                               kernel_type kernel) const {
  const LibintSystem& system = *m_system_ptr;

  // TODO libint2::Engine is not thread-safe, so for running this in
  // parallel, we need one engine object for each thread.
  //
  // TODO Parallelise here!
  //
  // TODO Code duplication with the ERICores

  const size_t max_nprim     = system.max_nprim();  // Maximum number of primitives
  const int max_l            = system.max_l();      // Maximum angular momentum
  const int derivative_order = 0;                   // Calculate no derivatives
  const real_type tolerance  = std::numeric_limits<real_type>::epsilon();
  libint2::Engine int_engine(libint2::Operator::coulomb, max_nprim, max_l,
                             derivative_order, tolerance);

  const LibintBasisShellData data{system};
  const krims::Range<size_t> a_shells = data.containing_shell_range(block[0]);
  const krims::Range<size_t> b_shells = data.containing_shell_range(block[1]);
  const krims::Range<size_t> c_shells = data.containing_shell_range(block[2]);
  const krims::Range<size_t> d_shells = data.containing_shell_range(block[3]);

  for (size_t sa : a_shells) {
    for (size_t sb : b_shells) {
      for (size_t sc : c_shells) {
        for (size_t sd : d_shells) {
          const auto& basis = system.basis();

          // Compute integrals (a b | c d), i.e. as in chemists/Mullikan notation
          // the shells a and b are basis functions of the same electron and so
          // are c and d.
          const size_t derivative_order = 0;
          int_engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                              derivative_order>(basis[sa], basis[sb], basis[sc],
                                                basis[sd]);

          // Extract computed values. If all values are zero (nullptr returned)
          // Skip this shell pair set alltogether.
          const scalar_type* values = int_engine.results()[0];
          if (values == nullptr) continue;

          // Call the kernel to do further computation
          kernel({{data.range_bfct(sa), data.range_bfct(sb), data.range_bfct(sc),
                   data.range_bfct(sd)}},
                 values);
        }  // sd
      }    // sc
    }      // sb
  }        // sa
}

//
// LibintCollection
//
const std::string IntegralCollection::id = "gaussian/libint";

Basis construct_basis(const krims::GenMap& parameters) {
  // TODO Generalise this: This function is pretty much the same for
  //      *all* gaussian basis set libraries.

  // The key: Either a file or the name of a basis set
  const std::string& key =
        parameters.at<const std::string>(IntegralLookupKeys::basis_set_name);

  const std::string errmsg =
        " Could not construct gaussian basis for the structure (provided via key \"" +
        IntegralLookupKeys::structure + "\") and the gaussian basis set \"" + key +
        "\" (provided via key \"" + IntegralLookupKeys::basis_set_name + "\")." +
        "\n\nDetails:\n--------\n";

  // If key is a valid file => Read it directly,
  // else interpret it as a basis set name and look it up in the library
  // of basis set files shipped with gint.
  const BasisSet basis_set = [&key, &errmsg] {
    try {
      return krims::path_exists(key) ? read_basisset(key) : lookup_basisset(key);
    } catch (const krims::ExcDatafileNotFound& e) {
      assert_throw(false, ExcInvalidIntegralParameters(errmsg + e.extra()));
    } catch (const ExcInvalidBasisSetFile& e) {
      assert_throw(false, ExcInvalidIntegralParameters(errmsg + e.extra()));
    }
    return BasisSet();
  }();

  auto structure_ptr = parameters.at_ptr<const Structure>(IntegralLookupKeys::structure);
  try {
    return Basis{*structure_ptr, basis_set};
  } catch (const ExcNoBasisForAtom& e) {
    assert_throw(false, ExcInvalidIntegralParameters(errmsg + e.extra()));
    return Basis{};
  }
}

LibintSystem make_system(const krims::GenMap& parameters) {
  auto structure_ptr = parameters.at_ptr<const Structure>(IntegralLookupKeys::structure);

  // Check if an explicit basis is specified: Use it
  // else construct a basis from the parameters
  if (parameters.exists(IntegralLookupKeys::basis)) {
    return LibintSystem(structure_ptr,
                        parameters.at<const Basis>(IntegralLookupKeys::basis));
  } else {
    return LibintSystem(structure_ptr, construct_basis(parameters));
  }
}

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : m_system{make_system(parameters)}, m_global{}, m_eri_tensor{m_system, m_global} {
  assert_throw(
        m_system.max_l() <= LIBINT_MAX_AM,
        ExcInvalidIntegralParameters(
              "Angular momentum limit exceeded for using libint. Angular momentum "
              "requested == " +
              std::to_string(m_system.max_l()) +
              " vs. max. angular momentum supported == " + std::to_string(LIBINT_MAX_AM) +
              ". Choose a smaller basis set or recompilie libint with a larger angular "
              "momentum limit."));
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  using libint2::Operator;

  switch (type) {
    case IntegralType::nuclear_attraction:
      return make_integral<OneElecIntegralCore>(Operator::nuclear, m_system, m_global);
    case IntegralType::overlap:
      return make_integral<OneElecIntegralCore>(Operator::overlap, m_system, m_global);
    case IntegralType::kinetic:
      return make_integral<OneElecIntegralCore>(Operator::kinetic, m_system, m_global);

    case IntegralType::coulomb: /* nobreak */
    case IntegralType::exchange:
      return make_integral<ERICore>(type, m_system, m_global);
  }

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<stored_matrix_type>(nullptr);
}

//
// LibintIntegralCoreBase
//

scalar_type LibintIntegralCoreBase::operator()(size_t row, size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  scalar_type ret = 0;
  auto kernel     = [&row, &col, &ret](LibintShell s, LibintShell t,
                                   const scalar_type* values) {
    assert_internal(s.first_bfct <= row);
    assert_internal(t.first_bfct <= col);

    if (values == nullptr) return;  // All values are zero => ret will be zero, too

    // Find the one value we need and extract it.
    const size_t fs = row - s.first_bfct;
    const size_t ft = col - t.first_bfct;
    ret             = values[fs * t.n_bfct + ft];
  };

  // Compute the one shell pair of the kernel we truly need:
  compute({row, row + 1}, {col, col + 1}, std::move(kernel));
  return ret;
}

void LibintIntegralCoreBase::extract_block(stored_matrix_type& M, const size_t start_row,
                                           const size_t start_col,
                                           const lazyten::Transposed mode,
                                           const scalar_type c_A,
                                           const scalar_type c_M) const {
  assert_finite(c_A);
  assert_finite(c_M);
  assert_greater_equal(start_row + M.n_rows(), n_rows());
  assert_greater_equal(start_col + M.n_cols(), n_cols());
  assert_sufficiently_tested(mode != lazyten::Transposed::ConjTrans);

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
  const auto col_range = krims::range(start_col, start_col + M.n_cols());

  auto kernel = [&row_range, &col_range, &M, &c_A](LibintShell s, LibintShell t,
                                                   const scalar_type* values) {
    if (values == nullptr) return;  // The whole shell-pair was screened away

    // The first basis function indices we care about in this shell pair block:
    size_t i_beg = std::max(s.first_bfct, row_range.lower_bound());
    size_t j_beg = std::max(t.first_bfct, col_range.lower_bound());

    // The last basis function indices we a care about in this shell pair block
    const size_t i_end = std::min(s.first_bfct + s.n_bfct, row_range.upper_bound());
    const size_t j_end = std::min(t.first_bfct + t.n_bfct, col_range.upper_bound());

    // libint2 packs the integrals into the shellsets in row-major
    // (C-like) form, so this iterates over the results and copy to the matrix M
    for (size_t i = i_beg; i < i_end; ++i) {
      for (size_t j = j_beg; j < j_end; ++j) {
        // Compute offset into values array:
        const size_t fs  = i - s.first_bfct;
        const size_t ft  = j - t.first_bfct;
        const size_t fst = fs * t.n_bfct + ft;

        // Shift i and j for placement into M:
        const size_t i_shifted = i - row_range.lower_bound();
        const size_t j_shifted = j - col_range.lower_bound();
        assert_internal(fst < s.n_bfct * t.n_bfct);

        M(i_shifted, j_shifted) += c_A * values[fst];
      }  // j
    }    // i
  };
  compute(row_range, col_range, std::move(kernel));
}

void LibintIntegralCoreBase::apply(const const_multivector_type& x, multivector_type& y,
                                   const lazyten::Transposed mode, const scalar_type c_A,
                                   const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != lazyten::Transposed::ConjTrans);
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

  auto kernel = [&x, &y, &c_A](LibintShell s, LibintShell t, const scalar_type* values) {
    if (values == nullptr) return;  // The whole shell-pair was screened away

    // libint2 packs the integrals into the shellsets in row-major
    // (C-like) form, so this iterates over the results one after the other
    // and applies it to the vectors
    for (size_t fs = 0; fs < s.n_bfct; ++fs) {
      for (size_t ft = 0; ft < t.n_bfct; ++ft) {
        const size_t fst = fs * t.n_bfct + ft;  // Index of function pair in result
        const size_t i   = s.first_bfct + fs;   // Absolute index in first shell
        const size_t j   = t.first_bfct + ft;   // Absolute index in second shell

        for (size_t vec = 0; vec < x.n_cols(); ++vec) {
          y(i, vec) += c_A * values[fst] * x(j, vec);
        }  // vec
      }    // j
    }      // i
  };
  compute(std::move(kernel));
}

//
// OneElecIntegralCore
//

void OneElecIntegralCore::compute(const krims::Range<size_t>& rows,
                                  const krims::Range<size_t>& cols,
                                  typename base_type::kernel_type&& kernel) const {
  const LibintBasisShellData data{base_type::system()};
  const krims::Range<size_t> row_shells = data.containing_shell_range(rows);
  const krims::Range<size_t> col_shells = data.containing_shell_range(cols);

  // TODO libint2::Engine is not thread-safe, so for running this in
  // parallel, we need one engine object for each thread.
  //
  // TODO Parallelise here!

  const LibintSystem& system = base_type::system();
  const size_t max_nprim     = system.max_nprim();  // Maximum number of primitives
  const int max_l            = system.max_l();      // Maximum angular momentum
  const int derivative_order = 0;                   // Calculate no derivatives
  const real_type tolerance  = std::numeric_limits<real_type>::epsilon();
  libint2::Engine int_engine(m_operator, max_nprim, max_l, derivative_order, tolerance);

  // If this is a nuclear attraction operator, set the point charges:
  if (m_operator == libint2::Operator::nuclear) {
    int_engine.set_params(base_type::system().point_charges());
  }

  // Loop over shell sets {s,t}
  for (size_t sa : row_shells) {
    for (size_t sb : col_shells) {
      // Compute the integral for a shell pair and call the kernel function with the
      // location where the computed integrals are stored.
      const auto& basis = system.basis();
      int_engine.compute1(basis[sa], basis[sb]);
      kernel(data.shell_info(sa), data.shell_info(sb), int_engine.results()[0]);
    }  // sb
  }    // sa
}

//
// ERI Integral core
//

void ERICore::compute(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
                      typename base_type::kernel_type&& kernel) const {
  const LibintBasisShellData data{base_type::system()};
  const krims::Range<size_t> row_shells = data.containing_shell_range(rows);
  const krims::Range<size_t> col_shells = data.containing_shell_range(cols);

  // TODO libint2::Engine is not thread-safe, so for running this in
  // parallel, we need one engine object for each thread.
  //
  // TODO Parallelise here!
  //
  // TODO Code duplication with the ERITensor

  const LibintSystem& system = base_type::system();
  const size_t max_nprim     = system.max_nprim();  // Maximum number of primitives
  const int max_l            = system.max_l();      // Maximum angular momentum
  const int derivative_order = 0;                   // Calculate no derivatives
  const real_type tolerance  = std::numeric_limits<real_type>::epsilon();
  libint2::Engine int_engine(libint2::Operator::coulomb, max_nprim, max_l,
                             derivative_order, tolerance);

  // Compute density matrix:
  const auto dens = compute_density_matrix();

  // Loop over shell sets {sa,sb,sc,sd}
  for (size_t sa : row_shells) {
    for (size_t sb : col_shells) {
      // Allocate a buffer to contain the values in the shell pair {sa,sb}
      // (stored in row-major)
      std::vector<scalar_type> buffer(data.n_bfct(sa) * data.n_bfct(sb), 0.);

      // Flag to indicate that all computed values of the inner double loop are zero
      bool all_zero = true;

      for (size_t sc = 0; sc < data.n_shells(); ++sc) {
        for (size_t sd = 0; sd < data.n_shells(); ++sd) {
          // Swap shell indices sa and sc if computing exchange
          const bool exchange = m_type == IntegralType::exchange;
          size_t sA           = exchange ? sc : sa;
          size_t sC           = exchange ? sa : sc;
          const auto& basis   = system.basis();

          // Compute integrals (A b | C d), i.e. as in chemists/Mullikan notation
          // the shells A and b are on the same centre, so are C and d.
          const size_t derivative_order = 0;
          int_engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                              derivative_order>(basis[sA], basis[sb], basis[sC],
                                                basis[sd]);

          // Extract computed values. If all values are zero (nullptr returned)
          // Skip this shell pair set alltogether.
          const scalar_type* values = int_engine.results()[0];
          if (values == nullptr) continue;

          all_zero = false;
          for (size_t a = 0; a < data.n_bfct(sa); ++a) {
            for (size_t b = 0; b < data.n_bfct(sb); ++b) {
              const size_t i_ab = a * data.n_bfct(sb) + b;

              for (size_t c = 0; c < data.n_bfct(sc); ++c) {
                // Swap shell indices a and c for exchange
                size_t A = exchange ? c : a;
                size_t C = exchange ? a : c;

                // Compute partial offsets
                size_t i_Ab  = A * data.n_bfct(sb) + b;
                size_t i_AbC = i_Ab * data.n_bfct(sC) + C;

                for (size_t d = 0; d < data.n_bfct(sd); ++d) {
                  // Indices to access element dens(c,d) of the density matrix
                  const size_t idens_c = data.first_bfct(sc) + c;
                  const size_t idens_d = data.first_bfct(sd) + d;

                  // Compute index into integral values:
                  const size_t i_AbCd = i_AbC * data.n_bfct(sd) + d;

                  // Perform explicit bound checks:
                  assert_internal(i_ab < buffer.size());
                  assert_internal(i_AbCd < data.n_bfct(sa) * data.n_bfct(sb) *
                                                 data.n_bfct(sc) * data.n_bfct(sd));

                  buffer[i_ab] += dens(idens_c, idens_d) * values[i_AbCd];
                }  // d
              }    // c
            }      // b
          }        // a

        }  // sh_d
      }    // sh_c

      if (all_zero) {
        kernel(data.shell_info(sa), data.shell_info(sb), nullptr);
      } else {
        kernel(data.shell_info(sa), data.shell_info(sb), buffer.data());
      }
    }  // sh_b
  }    // sh_a
}

}  // namespace libint
}  // namespace gaussian
}  // namespace gint
#endif  // GINT_HAVE_LIBINT
