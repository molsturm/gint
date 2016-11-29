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

void apply_stored_matrix(const real_stored_mtx_type& A,
			const_multivector_type & x,
			multivector_type& y,
			const linalgwrap::Transposed mode,
			const scalar_type c_A, const scalar_type c_y)  {
  for (size_t i = 0; i < x.n_rows(); i++)
    for (size_t j = 0; j < x.n_cols(); j++)
      y(i,j) = (c_y == 0? 0 : c_y * y(i,j));    

  // Everything in this module is real and symmetric, so we can ignore mode.
  for(size_t i=0;i<A.n_rows();i++)
    for(size_t j=0;j<A.n_cols();j++)
      for(size_t k=0;k<x.n_cols();k++)
	y(i,k) += c_A * A(i,j) * x(i,k);
}
  
  
}  // namespace static14
}  // namespace atomic
}  // namespace gint
  

#endif // ifdef GINT_STATIC_INTEGRALS
