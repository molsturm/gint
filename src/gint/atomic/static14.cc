#ifdef GINT_STATIC_INTEGRALS
#include "static14.hh"
//
// Implementation
//
namespace gint {
namespace atomic {
namespace static14 {

  
const std::string IntegralCollection::id = "atomic/static14";
const std::string IntegralCollection::name = "Fully precomputed integral data for 14 sturmians";


IntegralCollection::IntegralCollection(const krims::ParameterMap& parameters)
  : k_exponent{parameters.at<double>("k_exponent")},
    Z_charge{parameters.at<double>("Z_charge")} {}

Integral<real_stored_mtx_type>
IntegralCollection::lookup_integral(const std::string& integral_name) const {
  // TODO provide these keys as static constants somewhere
  if (integral_name == "nuclear_attraction")
    return make_integral<NuclearAttractionIntegralCore>(
          k_exponent, Z_charge);
  if (integral_name == "overlap")
    return make_integral<OverlapIntegralCore>();
  if (integral_name == "kinetic")
    return make_integral<KineticIntegralCore>(k_exponent);
  if (integral_name == "coulomb")
    return make_integral<ERICore>(false, k_exponent);
  if (integral_name == "exchange")
    return make_integral<ERICore>(true, k_exponent);

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<real_stored_mtx_type>(nullptr);
}

}  // namespace static14
}  // namespace atomic
}  // namespace gint
  

#endif // ifdef GINT_STATIC_INTEGRALS
