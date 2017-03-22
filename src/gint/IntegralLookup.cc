#include "IntegralLookup.hh"
#include <sstream>

#include "atomic/cs_dummy.hh"
#include "atomic/cs_naive.hh"
#include "atomic/cs_reference.hh"
#include "atomic/cs_reference_pc.hh"
#include "atomic/cs_static14.hh"
#include "gaussian/libint.hh"

namespace gint {

/** Function which registers the default basis types, which are present
 *  in this implementation of gint by default. */
void register_gint_basis_types() {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedef"
  //
  // Real atomic
  //
  typedef IntegralLookup<OrbitalType::REAL_ATOMIC> ra_t;

  //
  // Complex atomic
  //
  typedef IntegralLookup<OrbitalType::COMPLEX_ATOMIC> ca_t;
#ifdef GINT_HAVE_STATIC_INTEGRALS
  ca_t::register_basis_type(atomic::cs_static14::IntegralCollection::id,
                            atomic::cs_static14::IntegralCollection::create);
#endif  // GINT_HAVE_STATIC_INTEGRALS
  ca_t::register_basis_type(atomic::cs_dummy::IntegralCollection::id,
                            atomic::cs_dummy::IntegralCollection::create);
  ca_t::register_basis_type(atomic::cs_naive::IntegralCollection::id,
                            atomic::cs_naive::IntegralCollection::create);
  ca_t::register_basis_type(atomic::cs_reference::IntegralCollection::id,
                            atomic::cs_reference::IntegralCollection::create);
  ca_t::register_basis_type(atomic::cs_reference_pc::IntegralCollection::id,
                            atomic::cs_reference_pc::IntegralCollection::create);

  //
  // Real molecular
  //
  typedef IntegralLookup<OrbitalType::REAL_MOLECULAR> rm_t;
#ifdef GINT_HAVE_LIBINT
  rm_t::register_basis_type(gaussian::libint::IntegralCollection::id,
                            gaussian::libint::IntegralCollection::create);
#endif  // GINT_HAVE_LIBINT

  //
  // Complex molecular
  //
  typedef IntegralLookup<OrbitalType::COMPLEX_MOLECULAR> cm_t;
#pragma GCC diagnostic pop
}

/** Once flag, which makes sure that the registration function above is only
 *  called once */
std::once_flag once_register_gint_basis_types;

template <OrbitalType otype>
std::map<std::string, collection_generator_type<otype>>
      IntegralLookup<otype>::map_basis_collection_generator{};

template <OrbitalType otype>
IntegralLookup<otype>::IntegralLookup(const krims::GenMap& parameters) {
  std::call_once(once_register_gint_basis_types, register_gint_basis_types);

  std::string basis_type =
        parameters.at("basis_type", std::string("<No basis_type supplied>"));

  auto it = map_basis_collection_generator.find(basis_type);

  if (it == std::end(map_basis_collection_generator)) {
    std::stringstream list;
    for (const auto& kv : map_basis_collection_generator) {
      list << "'" << kv.first << "' ";
    }
    assert_throw(false,
                 ExcInvalidIntegralParameters(
                       "There is no integral library corresponding to basis_type='" +
                       basis_type + "'. Only " + list.str() + " are known."));
  }

  m_integral_collection_ptr = it->second(parameters);
}

// Explicitly instantiate the parts we need:
template class IntegralLookup<OrbitalType::COMPLEX_ATOMIC>;
template class IntegralLookup<OrbitalType::REAL_ATOMIC>;
template class IntegralLookup<OrbitalType::REAL_MOLECULAR>;
// TODO Does not compile at the moment:
//   template class IntegralLookup<OrbitalType::COMPLEX_MOLECULAR>;

}  // namespace gint
