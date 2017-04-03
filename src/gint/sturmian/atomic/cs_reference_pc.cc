#include "cs_reference_pc.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_reference_pc {

const std::string IntegralCollection::id = "atomic/cs_reference_pc";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : m_system{}, m_integral_calculator{} {
  if (parameters.exists("nlmbasis")) {
    m_system.basis = parameters.at<const NlmBasis>("nlmbasis");
  } else {
    const int n_max = parameters.at<int>("n_max");
    const int l_max = parameters.at<int>("l_max");
    const int m_max = parameters.at<int>("m_max");

    m_system.basis = NlmBasis(n_max, l_max, m_max);
  }

  m_system.Z = parameters.at<scalar_type>("Z_charge");
  m_system.k = parameters.at<scalar_type>("k_exponent");
  m_integral_calculator = sturmint::atomic::cs_reference_pc::Atomic(m_system.basis);
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

scalar_type ERICore::operator()(size_t a, size_t b) const {
  assert_greater(a, n_rows());
  assert_greater(b, n_cols());
  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  // Compute density matrix:
  const auto density = linalgwrap::outer_prod_sum(coeff_bo(), coeff_bo());

  real_type sum = 0;
  size_t norb = m_integral_calculator.n_bas();
  for (size_t c = 0; c < norb; c++) {
    // Swap b and c if computing exchange
    const bool exchange = type() == IntegralType::exchange;
    size_t A = exchange ? c : a;
    size_t C = exchange ? a : c;

    for (size_t d = 0; d < norb; d++)
      sum += m_integral_calculator.repulsion(A, b, C, d) * density(c, d);
  }
  return system().k * sum;
}

}  // namespace cs_reference_pc
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
