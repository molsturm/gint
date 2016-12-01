#include "IntegralLookup.hh"
#include <sstream>

namespace gint {

template <>
IntegralLookup<COMPLEX_ATOMIC>::IntegralLookup(const krims::ParameterMap& parameters) {
  std::string basis_type = parameters.at("basis_type", std::string("No basis_type supplied"));

  const auto& create_collection(basis_type_map_ca.find(basis_type));

#if DEBUG
  std::stringstream list;
  for (const auto& kv : basis_type_map_ca) {
    list << "'" << kv.first << "' ";
  }
  assert_throw(
        create_collection != basis_type_map_ca.end(),
        ExcInvalidIntegralParameters("There is no integral library corresponding to basis_type='" +
                                     basis_type + "'. Only " + list.str() + " are known."));
#endif

  m_integral_collection = create_collection->second(parameters);
}

}  // namespace gint
