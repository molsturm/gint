#include "IntegralLookup.hh"
#include "OrbitalType.hh"
#include <mutex>
#include <sstream>

#include "atomic/cs_dummy.hh"
#include "atomic/cs_naive.hh"
#include "atomic/cs_reference.hh"
#include "atomic/cs_reference_pc.hh"
#include "atomic/cs_static14.hh"
#include "gaussian/libint.hh"

namespace gint {

namespace {
#define REGISTER_BASIS_TYPE(LOOKUP, BASIS)                          \
  {                                                                 \
    LOOKUP::register_basis_type(BASIS::IntegralCollection::id,      \
                                BASIS::IntegralCollection::create); \
  }

/** Function which registers the default basis types, which are present
 *  in this implementation of gint by default. */
void register_gint_basis_types() {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
  typedef IntegralLookup<real_valued::stored_matrix_type> real_t;
  typedef IntegralLookup<complex_valued::stored_matrix_type> complex_t;

// Coulomb-Sturmians
#ifdef GINT_HAVE_STATIC_INTEGRALS
  REGISTER_BASIS_TYPE(real_t, atomic::cs_static14);
#endif  // GINT_HAVE_STATIC_INTEGRALS
  REGISTER_BASIS_TYPE(real_t, atomic::cs_dummy);
  REGISTER_BASIS_TYPE(real_t, atomic::cs_naive);
  REGISTER_BASIS_TYPE(real_t, atomic::cs_reference);
  REGISTER_BASIS_TYPE(real_t, atomic::cs_reference_pc);

// Gaussians
#ifdef GINT_HAVE_LIBINT
  REGISTER_BASIS_TYPE(real_t, gaussian::libint);
#endif  // GINT_HAVE_LIBINT

#pragma GCC diagnostic pop  // -Wunused-local-typedef
}
#undef REGISTER_BASIS_TYPE

/** Once flag, which makes sure that the registration function above is only
 *  called once */
std::once_flag once_register_gint_basis_types;
}  // namespace

template <typename StoredMatrix>
std::map<std::string, collection_generator_type<StoredMatrix>>
      IntegralLookup<StoredMatrix>::map_basis_collection_generator{};

template <typename StoredMatrix>
IntegralLookup<StoredMatrix>::IntegralLookup(const krims::GenMap& parameters) {
  std::call_once(once_register_gint_basis_types, register_gint_basis_types);

  // TODO The orbital type is entirely ignored at the moment.
  //      We would need different integral collections for the different orbital types.
  //      Think about a way to build this into the existing system when we need it.
  const OrbitalType otype =
        parameters.at<OrbitalType>("orbital_type", OrbitalType::REAL_MOLECULAR);
  (void)otype;

  const auto basis_type =
        parameters.at<std::string>("basis_type", "<No basis_type supplied>");
  auto it = map_basis_collection_generator.find(basis_type);

  if (it == std::end(map_basis_collection_generator)) {
    std::stringstream list;
    for (const auto& kv : map_basis_collection_generator) {
      list << "'" << kv.first << "' ";
    }
    assert_throw(false,
                 ExcInvalidIntegralParameters(
                       "There is no integral collection corresponding to basis_type='" +
                       basis_type + "'. Only " + list.str() + " are known."));
  }

  m_integral_collection_ptr = it->second(parameters);
}

// Explicitly instantiate real and complex version:
template class IntegralLookup<real_valued::stored_matrix_type>;
template class IntegralLookup<complex_valued::stored_matrix_type>;

}  // namespace gint
