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
      : IntegralCollectionBase(parameters), m_integral_calculator(m_system.basis) {}

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
      assert_implemented(false);
      return Integral<stored_matrix_type>(nullptr);
  }
}

}  // namespace cs_reference
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
