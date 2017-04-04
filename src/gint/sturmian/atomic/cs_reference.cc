#include "cs_reference.hh"
#include "nlm_order/ERICore.hh"
#include "nlm_order/ERICoreHighPrecision.hh"
#include "nlm_order/OneElectronIntegralCores.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_reference {

const std::string IntegralCollection::id = "atomic/cs_reference";

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
  m_integral_calculator = sturmint::atomic::cs_reference::Atomic(m_system.basis);
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  typedef ERICore<Atomic> eri_type;
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

}  // namespace cs_reference
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
