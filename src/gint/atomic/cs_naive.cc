#include "cs_naive.hh"

namespace gint {
namespace atomic {
namespace cs_naive {

const std::string IntegralCollection::id = "atomic/cs_naive";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{parameters.at<double>("Z_charge")},
        basis{parameters.at<int>("n_max"), parameters.at<int>("l_max"),
              parameters.at<int>("m_max")},
        integral_calculator{basis} {
  int n_max = parameters.at<int>("n_max");
  int l_max = parameters.at<int>("l_max");
  int m_max = parameters.at<int>("m_max");

  assert_throw(0 < n_max && n_max <= 6,
               ExcInvalidIntegralParameters(
                     "Maximum principle quantum number (" + std::to_string(n_max) +
                     ") needs to be in the range [1,6] for cs_naive, since higher "
                     "values are not yet implemented."));
  assert_throw(0 <= l_max && l_max <= 4,
               ExcInvalidIntegralParameters(
                     "Maximum angular momentum quantum number (" + std::to_string(l_max) +
                     ") needs to be in the range [0,4] for cs_naive, since higher values "
                     "are not yet implemented."));
  assert_throw(
        0 <= m_max && m_max <= 4,
        ExcInvalidIntegralParameters("Maximum magnetic momentum quantum number (" +
                                     std::to_string(m_max) +
                                     ") needs to be in the range [0,4] for cs_naive, "
                                     "since higher values are not yet implemented."));
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  switch (type) {
    case IntegralType::nuclear_attraction:
      return make_integral<NuclearAttractionIntegralCore>(integral_calculator, k_exponent,
                                                          Z_charge);
    case IntegralType::overlap:
      return make_integral<OverlapIntegralCore>(integral_calculator);
    case IntegralType::kinetic:
      return make_integral<KineticIntegralCore>(integral_calculator, k_exponent);
    case IntegralType::coulomb:
      return make_integral<ERICore>(integral_calculator, false, k_exponent);
    case IntegralType::exchange:
      return make_integral<ERICore>(integral_calculator, true, k_exponent);
  }

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<stored_matrix_type>(nullptr);
}

//
// ERI Integral core
//

void ERICore::apply(const const_multivector_type& x, multivector_type& y,
                    const linalgwrap::Transposed mode, const scalar_type c_A,
                    const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  const size_t norb = x.n_rows(), n_vectors = x.n_cols();

  for (size_t i = 0; i < y.n_rows(); i++)
    for (size_t j = 0; j < y.n_cols(); j++) y(i, j) = (c_y != 0 ? c_y * y(i, j) : 0);

  // Work internally in highest available precision
  vector<data_real_t> ysum(n_vectors);
  for (size_t a = 0; a < norb; a++) {
    memset(&ysum[0], 0, n_vectors * sizeof(sturmint::data_real_t));
    for (size_t b = 0; b < norb; b++) {
      data_real_t JKab = (*this)(a, b);

      for (size_t q = 0; q < n_vectors; q++) {
        ysum[q] += c_A * JKab * x(b, q);
      }
    }
    for (size_t q = 0; q < n_vectors; q++) y(a, q) += ysum[q];
  }
}

scalar_type ERICore::operator()(size_t a, size_t b) const {
  assert_greater(a, n_rows());
  assert_greater(b, n_cols());
  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  const coefficients_type& Cocc(*coefficients_occupied_ptr);
  stored_matrix_type density(basis.size(), basis.size());
  for (size_t p = 0; p < coefficients_occupied_ptr->n_vectors(); p++) {
    const auto& C = Cocc[p];
    for (size_t c = 0; c < basis.size(); c++)
      for (size_t d = 0; d < basis.size(); d++) density(c, d) += C[c] * C[d];
  }

  data_real_t sum = 0;

  for (size_t c = 0; c < basis.size(); c++) {
    size_t A = exchange ? c : a;  // Swap b and c if computing exchange // TODO: Check
    size_t C = exchange ? a : c;

    for (size_t d = 0; d < basis.size(); d++) {
      sturmint::data_real_t density_cd = 0;
      for (size_t p = 0; p < Cocc.n_vectors(); p++) density_cd += Cocc[p][c] * Cocc[p][d];

      sum += m_integral_calculator.repulsion(basis[A], basis[b], basis[C], basis[d]) *
             density_cd;
    }
  }
  return static_cast<scalar_type>(k * sum);
}

void ERICore::update(const krims::GenMap& map) {
  const std::string occ_coeff_key = Integral<stored_matrix_type>::update_key_coefficients;

  if (!map.exists(occ_coeff_key)) return;

  // Get coefficients as a shared pointer (having ownership)
  coefficients_occupied_ptr =
        static_cast<coefficients_ptr_type>(map.at_ptr<coefficients_type>(occ_coeff_key));

  // We will contract the coefficient row index over the number of
  // basis functions.
  if (coefficients_occupied_ptr->n_vectors() == 0) return;
}

}  // namespace cs_naive
}  // namespace atomic
}  // namespace gint
