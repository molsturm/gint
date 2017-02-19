#include "cs_dummy.hh"

namespace gint {
namespace atomic {
namespace cs_dummy {

// ----------------------------------------------------------------------
//			    IMPLEMENTATION
// ----------------------------------------------------------------------
std::string IntegralCollection::id = "atomic/cs_dummy";
std::string IntegralCollection::name = "Dummy implementation of atomic Coulomb Sturmians";

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
      const std::string& integral_name) const {
  // TODO provide these keys as static constants somewhere
  if (integral_name == "nuclear_attraction")
    return make_integral<NuclearAttractionIntegralCore>(integral_calculator, k_exponent,
                                                        Z_charge);
  if (integral_name == "overlap")
    return make_integral<OverlapIntegralCore>(integral_calculator);
  if (integral_name == "kinetic")
    return make_integral<KineticIntegralCore>(integral_calculator, k_exponent);
  if (integral_name == "coulomb")
    return make_integral<ERICore>(integral_calculator, false, k_exponent);
  if (integral_name == "exchange")
    return make_integral<ERICore>(integral_calculator, true, k_exponent);

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<real_stored_mtx_type>(nullptr);
}

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace gint
