#include "cs_naive.hh"
#include "nlm_order/ERICore.hh"
#include "nlm_order/ERICoreHighPrecision.hh"
#include "nlm_order/OneElectronIntegralCores.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_naive {

const std::string IntegralCollection::id = "atomic/cs_naive";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : m_system{}, m_integral_calculator{} {
  const int n_max = parameters.at<int>("n_max");
  const int l_max = parameters.at<int>("l_max");
  const int m_max = parameters.at<int>("m_max");

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

  m_system.Z = parameters.at<scalar_type>("Z_charge");
  m_system.k = parameters.at<scalar_type>("k_exponent");
  m_system.basis = NlmBasis(n_max, l_max, m_max);
  m_integral_calculator = Atomic(m_system.basis);
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  typedef ERICoreHighPrecision<Atomic> eri_type;
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
      return make_integral<eri_type>(m_integral_calculator, m_system,
                                     IntegralIdentifier{id, type});

    default:
      assert_dbg(false, krims::ExcNotImplemented());
      return Integral<stored_matrix_type>(nullptr);
  }
}

}  // namespace cs_naive
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
