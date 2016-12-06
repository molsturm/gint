#include "cs_naive.hh"


namespace gint {
namespace atomic {
namespace cs_naive {

// ----------------------------------------------------------------------
//			    IMPLEMENTATION
// ----------------------------------------------------------------------
std::string IntegralCollection::id   = "cs_atomic/naive";
std::string IntegralCollection::name = "Naive implementation of atomic Coulomb Sturmians";

IntegralCollection::IntegralCollection(
      const krims::ParameterMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{parameters.at<double>("Z_charge")},
	basis{parameters.at<int>("n_max"),parameters.at<int>("l_max"),parameters.at<int>("m_max")},
        integral_calculator{basis} {
	  //	  cerr << "Built integral collection, number of orbitals is " << basis.length << "\n";
	}

Integral<real_stored_mtx_type> IntegralCollection::lookup_integral(
      const std::string& integral_name) const {
  // TODO provide these keys as static constants somewhere
  if (integral_name == "nuclear_attraction")
    return make_integral<NuclearAttractionIntegralCore>(integral_calculator, k_exponent, Z_charge);
  if (integral_name == "overlap")  return make_integral<OverlapIntegralCore>(integral_calculator);
  if (integral_name == "kinetic")  return make_integral<KineticIntegralCore>(integral_calculator, k_exponent);
  if (integral_name == "coulomb")  return make_integral<ERICore>(integral_calculator, false, k_exponent);
  if (integral_name == "exchange") return make_integral<ERICore>(integral_calculator, true, k_exponent);

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<real_stored_mtx_type>(nullptr);
}

} // namespace cs_naive
} // namespace atomic
} // namespace gint
