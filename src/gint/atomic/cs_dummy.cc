#include "cs_dummy.hh"

namespace gint {
namespace atomic {
namespace cs_dummy {

const std::string IntegralCollection::id = "atomic/cs_dummy";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{parameters.at<double>("Z_charge")},
        n_max{parameters.at<int>("n_max")},
        integral_calculator{n_max} {
  assert_throw(n_max > 0, ExcInvalidIntegralParameters(
                                "Maximum principle quantum number (" +
                                std::to_string(n_max) + ") needs to be greater 0."));
  assert_throw(n_max <= 3,
               ExcInvalidIntegralParameters(
                     "cs_dummy is only implemented up to n_max==3. You provided "
                     "a maximum principle quantum number of " +
                     std::to_string(n_max) + ", which is too large."));
}

Integral<real_stored_mtx_type> IntegralCollection::lookup_integral(
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
  return Integral<real_stored_mtx_type>(nullptr);
}

//
// ERI Integral core
//

void ERICore::apply(const const_multivector_type& x, multivector_type& y,
                    const linalgwrap::Transposed mode, const scalar_type alpha,
                    const scalar_type beta) const {
  using namespace sturmint::orbital_index;
  using namespace sturmint::atomic::cs_dummy;

  assert_finite(alpha);
  assert_finite(beta);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

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
  using namespace sturmint::atomic::cs_dummy;
  using namespace sturmint::orbital_index;

  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  const coefficients_type& Cocc(*coefficients_occupied_ptr);

  stored_mtx_type density(norb, norb);
  for (size_t p = 0; p < coefficients_occupied_ptr->n_vectors(); p++) {
    const auto& C = Cocc[p];
    for (size_t c = 0; c < norb; c++)
      for (size_t d = 0; d < norb; d++) density(c, d) += C[c] * C[d];
  }

  real_type sum = 0;
  for (size_t c = 0; c < norb; c++) {
    size_t A = exchange ? c : a;  // Swap b and c if computing exchange // TODO: Check
    size_t C = exchange ? a : c;

    size_t i_abc = norb * (C + norb * (b + norb * A));

    for (size_t d = 0; d < norb; d++) sum += repulsion14[i_abc + d] * density(c, d);
  }
  return k * sum;
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
}

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace gint
