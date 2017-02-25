#ifdef GINT_HAVE_STATIC_INTEGRALS
#include "cs_static14.hh"

namespace gint {
namespace atomic {
namespace cs_static14 {

const std::string IntegralCollection::id = "atomic/cs_static14";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{parameters.at<double>("Z_charge")} {}

Integral<real_stored_mtx_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  switch (type) {
    case IntegralType::nuclear_attraction:
      return make_integral<NuclearAttractionIntegralCore>(k_exponent, Z_charge);
    case IntegralType::overlap:
      return make_integral<OverlapIntegralCore>();
    case IntegralType::kinetic:
      return make_integral<KineticIntegralCore>(k_exponent);
    case IntegralType::coulomb:
      return make_integral<ERICore>(false, k_exponent);
    case IntegralType::exchange:
      return make_integral<ERICore>(true, k_exponent);
  }

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<real_stored_mtx_type>(nullptr);
}

void apply_stored_matrix(const real_stored_mtx_type& A, const_multivector_type& x,
                         multivector_type& y, const linalgwrap::Transposed mode,
                         const scalar_type c_A, const scalar_type c_y) {
  // scale y by c_y or set to zero
  for (size_t i = 0; i < x.n_rows(); i++)
    for (size_t j = 0; j < x.n_cols(); j++) y(i, j) = (c_y == 0 ? 0 : c_y * y(i, j));

  // Everything in this module is real and symmetric, so we can ignore mode.
  for (size_t k = 0; k < x.n_cols(); k++)  // Iterate over vectors
    for (size_t i = 0; i < A.n_rows(); i++) {
      real_type row_sum = 0;
      for (size_t j = 0; j < A.n_cols(); j++) row_sum += c_A * A(i, j) * x(j, k);
      y(i, k) += row_sum;
    }
}

//
// IntegralCores
//

typename ERICore::scalar_type ERICore::operator()(size_t a, size_t b) const {
  assert_greater(a, n_rows());
  assert_greater(b, n_cols());

  // This routine computes
  //
  // J_{ab} = \sum_{cd} P_{cd} < ac | bd >
  //        = \sum_{cd} \sum_{k \in occ} (C^{k})_c (C^{k})_d < ac | bd >
  // where a and b are the same centre, so are c and d
  //
  // or alternatively
  //
  // K_{ab} = \sum_{cd} P_{cd} < ab | cd >
  // where a and c are the same centre, so are b and d

  // Repulsion integrals indexed in shell-pairs
  const auto& i_bbbb = detail::Static14Data<stored_mtx_type>::i_bbbb_base;
  const size_t nbas = detail::Static14Data<stored_mtx_type>::nbas;

  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  // Density matrix expression
  auto density_bb =
        outer_prod_sum(*coefficients_occupied_ptr, *coefficients_occupied_ptr);
  assert_dbg(density_bb.n_rows() == nbas, krims::ExcInternalError());
  assert_dbg(density_bb.n_cols() == nbas, krims::ExcInternalError());

  // Shell pair index for basis functions a and b:
  const size_t ab_pair = a * nbas + b;

  // Sum accumulator variable for this exchange or
  // coulomb  matrix element
  scalar_type mat_ab{0};

  // Double loop over basis functions c and d:
  for (size_t c = 0; c < nbas; ++c) {
    // Shell pair index for basis functions a and c:
    const size_t ac_pair = a * nbas + c;

    for (size_t d = 0; d < nbas; ++d) {
      // Shell pair index for basis functions c and d:
      // or b and d:
      const size_t cd_pair = c * nbas + d;
      const size_t db_pair = d * nbas + b;

      // Perform contraction:
      const scalar_type i_elem =
            exchange ? i_bbbb(ac_pair, db_pair) : i_bbbb(ab_pair, cd_pair);
      mat_ab += k * i_elem * density_bb(c, d);
    }  // d
  }    // c

  return mat_ab;
}

void ERICore::apply(const_multivector_type& x, multivector_type& y,
                    const linalgwrap::Transposed mode, const scalar_type c_A,
                    const scalar_type c_y) const {
  const size_t nbas = detail::Static14Data<stored_mtx_type>::nbas;
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), nbas);
  assert_size(y.n_rows(), nbas);
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  // Scale the current values of y or set them to zero
  // (if c_y == 0): We are now done with c_y and do not
  // need to worry about it any more in this function
  for (size_t i = 0; i < y.n_rows(); i++)
    for (size_t j = 0; j < y.n_cols(); j++) y(i, j) = (c_y != 0 ? c_y * y(i, j) : 0);

  // if c_this == 0 we are done
  if (c_A == linalgwrap::Constants<scalar_type>::zero) return;

  for (size_t veci = 0; veci < x.n_cols(); ++veci) {
    for (size_t row = 0; row < nbas; ++row) {
      scalar_type sum = 0;
      for (size_t k = 0; k < nbas; ++k) sum += (*this)(row, k) * x(k, veci);

      y(row, veci) += c_A * sum;
    }  // row
  }    // veci
}

void ERICore::extract_block(stored_matrix_type& M, const size_t start_row,
                            const size_t start_col, const linalgwrap::Transposed mode,
                            const scalar_type c_this, const scalar_type c_M) const {
  using namespace linalgwrap;
#ifdef DEBUG
  const size_t nbas = detail::Static14Data<stored_matrix_type>::nbas;
#endif

  assert_finite(c_this);
  assert_finite(c_M);
  assert_greater_equal(start_row + M.n_rows(), nbas);
  assert_greater_equal(start_col + M.n_cols(), nbas);
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  // Scale the current values of M or set them to zero
  // (if c_M == 0): We are now done with c_M and do not
  // need to worry about it any more in this function
  // TODO we are kind of calling an internal function here
  linalgwrap::detail::scale_or_set(M, c_M);

  if (c_this == Constants<scalar_type>::zero) return;

  for (size_t row = 0; row < M.n_rows(); ++row) {
    for (size_t col = 0; col < M.n_cols(); ++col) {
      switch (mode) {
        case Transposed::None:
        case Transposed::Trans:
          // No difference since symmetric
          M(row, col) += c_this * (*this)(start_row + row, start_col + col);
          break;
        case Transposed::ConjTrans:
          // A variant of std::conj, which does not return a complex
          // data type if scalar is real only.
          // TODO we are kind of calling an internal function here
          linalgwrap::detail::ConjFctr mconj;
          M(row, col) += c_this * mconj((*this)(start_col + col, start_row + row));
          break;
      }  // mode
    }    // col
  }      // row
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
  assert_size(coefficients_occupied_ptr->n_elem(),
              detail::Static14Data<stored_mtx_type>::nbas);
}

}  // namespace cs_static14
}  // namespace atomic
}  // namespace gint

#endif  // ifdef GINT_HAVE_STATIC_INTEGRALS
