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
    m_system.basis = parameters.at<const NlmBasis>("nlmbasis");
    m_repulsiondata_filename = parameters.at<string>("repulsiondata_filename");
  } else {
    const int nmax = parameters.at<int>("n_max");
    const int lmax = parameters.at<int>("l_max");
    const int mmax = parameters.at<int>("m_max");
    m_system.basis = NlmBasis(nmax, lmax, mmax);

    m_repulsiondata_filename = std::string("repulsiondata-nlm-") + to_string(nmax) + "-" +
                               to_string(lmax) + "-" + to_string(mmax) + ".bin";
  }

  m_system.Z = parameters.at<scalar_type>("Z_charge");
  m_system.k = parameters.at<scalar_type>("k_exponent");
  m_integral_calculator = Atomic(m_system.basis, m_repulsiondata_filename);
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

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
