#include "cs_dummy.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_dummy {

const std::string IntegralCollection::id = "atomic/cs_dummy";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : m_system{}, m_repulsiondata_filename{""}, m_integral_calculator{} {
  const bool explicit_basis =
        parameters.exists("nlmbasis") && parameters.exists("repulsiondata_filename");

  if (explicit_basis) {
    m_system.basis = parameters.at<const nlmCollection>("nlmbasis");
    m_repulsiondata_filename = parameters.at<string>("repulsiondata_filename");
  } else {
    const int nmax = parameters.at<int>("n_max");
    const int lmax = parameters.at<int>("l_max");
    const int mmax = parameters.at<int>("m_max");
    m_system.basis = nlmCollection(nmax, lmax, mmax);

    m_repulsiondata_filename = std::string("repulsiondata-nlm-") + to_string(nmax) + "-" +
                               to_string(lmax) + "-" + to_string(mmax) + ".bin";
  }

  m_system.Z = parameters.at<scalar_type>("Z_charge");
  m_system.k = parameters.at<scalar_type>("k_exponent");
  m_integral_calculator = Atomic(m_system.basis, m_repulsiondata_filename);
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  const std::string& id = IntegralCollection::id;

  switch (type) {
    case IntegralType::nuclear_attraction:
      return make_integral<NuclearAttractionIntegralCore>(m_system, id);
    case IntegralType::overlap:
      return make_integral<OverlapIntegralCore>(m_system, id);
    case IntegralType::kinetic:
      return make_integral<KineticIntegralCore>(m_system, id);

    case IntegralType::coulomb: /* nobreak */
    case IntegralType::exchange:
      return make_integral<ERICore>(m_integral_calculator, m_system, type);

    default:
      assert_dbg(false, krims::ExcNotImplemented());
      return Integral<stored_matrix_type>(nullptr);
  }
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
  assert_greater(a, n_rows());
  assert_greater(b, n_cols());
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
    // Swap b and c if computing exchange
    const bool exchange = type() == IntegralType::exchange;
    const size_t A = exchange ? c : a;
    const size_t C = exchange ? a : c;

    for (size_t d = 0; d < norb; d++)
      sum += m_integral_calculator.repulsion(A, b, C, d) * density(c, d);
  }
  return system().k * sum;
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

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
