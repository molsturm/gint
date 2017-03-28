#include "cs_reference.hh"

namespace gint {
namespace atomic {
namespace cs_reference {

const std::string IntegralCollection::id = "atomic/cs_reference";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{parameters.at<double>("Z_charge")} {
  if (parameters.exists("nlmbasis")) {
    basis = parameters.at<const vector<nlm_t>>("nlmbasis");
  } else {
    int nmax = parameters.at<int>("n_max");
    int lmax = parameters.at<int>("l_max");
    int mmax = parameters.at<int>("m_max");

    basis = nlmCollection(nmax, lmax, mmax);
  }

  integral_calculator = sturmint::atomic::cs_reference::Atomic(basis);
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
                    const linalgwrap::Transposed mode, const scalar_type alpha,
                    const scalar_type beta) const {
  using namespace sturmint::orbital_index;
  using namespace sturmint::atomic::cs_reference;

  assert_finite(alpha);
  assert_finite(beta);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  size_t norb = m_integral_calculator.n_bas();

  for (size_t i = 0; i < y.n_rows(); i++)
    for (size_t j = 0; j < y.n_cols(); j++) y(i, j) = (beta != 0 ? beta * y(i, j) : 0);

  for (size_t a = 0; a < norb; a++)
    for (size_t b = 0; b < norb; b++) {
      real_type JKab = (*this)(a, b);

      for (size_t q = 0; q < x.n_cols(); q++) {
        y(a, q) += alpha * JKab * x(b, q);
      }
    }
}

scalar_type ERICore::operator()(size_t a, size_t b) const {
  using namespace sturmint::atomic::cs_reference;
  using namespace sturmint::orbital_index;

  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  const coefficients_type& Cocc(*coefficients_occupied_ptr);

  size_t norb = m_integral_calculator.n_bas();
  stored_matrix_type density(norb, norb);
  for (size_t p = 0; p < coefficients_occupied_ptr->n_vectors(); p++) {
    const auto& C = Cocc[p];
    for (size_t c = 0; c < norb; c++)
      for (size_t d = 0; d < norb; d++) density(c, d) += C[c] * C[d];
  }

  real_type sum = 0;
  for (size_t c = 0; c < norb; c++) {
    size_t A = exchange ? c : a;  // Swap b and c if computing exchange // TODO: Check
    size_t C = exchange ? a : c;

    //    size_t i_abc = norb * (C + norb * (b + norb * A));

    //    for (size_t d = 0; d < norb; d++) sum +=
    //    m_integral_calculator.repulsionNxNxNxN[i_abc + d] * density(c, d);
    for (size_t d = 0; d < norb; d++)
      sum += m_integral_calculator.repulsion(A, b, C, d) * density(c, d);
  }
  return k * sum;
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
  assert_size(coefficients_occupied_ptr->n_elem(), m_integral_calculator.n_bas());
}

}  // namespace cs_reference
}  // namespace atomic
}  // namespace gint
