#include "IntegralLookup.hh"

namespace gint {
  
template<> 
IntegralLookup<COMPLEX_ATOMIC>::IntegralLookup(const krims::ParameterMap& parameters){
  std::string basis_type = parameters.at("basis_type", std::string("No basis_type supplied"));

  const auto &create_collection(basis_type_map_ca.find(basis_type));

  assert_throw(create_collection != basis_type_map_ca.end(),
               ExcInvalidIntegralParameters("There is no integral library corresponding to basis_type='"+basis_type+"'"));

  m_integral_collection = create_collection->second(parameters);
}

} // namespace gint