//
// Copyright (C) 2017 by the gint authors
//
// This file is part of gint.
//
// gint is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gint is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gint. If not, see <http://www.gnu.org/licenses/>.
//

#include "IntegralLookup.hh"
#include "IntegralLookupKeys.hh"
#include "OrbitalType.hh"
#include <mutex>
#include <sstream>

#include "gaussian/libcint.hh"
#include "gaussian/libint.hh"
#include "sturmian/atomic.hh"

namespace gint {

namespace {
#define REGISTER_BASIS_TYPE(LOOKUP, BASIS)                          \
  {                                                                 \
    LOOKUP::register_basis_type(BASIS::IntegralCollection::id,      \
                                BASIS::IntegralCollection::create); \
  }

typedef IntegralLookup<real_valued::stored_matrix_type> real_t;

/** Function which registers the default basis types, which are present
 *  in this implementation of gint by default. */
void register_gint_basis_types() {
// TODO: enable when needed:
//       typedef IntegralLookup<complex_valued::stored_matrix_type> complex_t;

// Coulomb-Sturmians
#ifdef GINT_HAVE_STATIC_INTEGRALS
  REGISTER_BASIS_TYPE(real_t, sturmian::atomic::cs_static14);
#endif  // GINT_HAVE_STATIC_INTEGRALS

#ifdef GINT_HAVE_STURMINT
  REGISTER_BASIS_TYPE(real_t, sturmian::atomic::cs_dummy);
  REGISTER_BASIS_TYPE(real_t, sturmian::atomic::cs_naive);
  REGISTER_BASIS_TYPE(real_t, sturmian::atomic::cs_reference);
  REGISTER_BASIS_TYPE(real_t, sturmian::atomic::cs_reference_pc);
#endif  // GINT_HAVE_STURMINT

// Gaussians
#ifdef GINT_HAVE_LIBINT
  REGISTER_BASIS_TYPE(real_t, gaussian::libint);
#endif  // GINT_HAVE_LIBINT
#ifdef GINT_HAVE_LIBCINT
  REGISTER_BASIS_TYPE(real_t, gaussian::libcint);
#endif  // GINT_HAVE_LIBCINT
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
std::vector<std::string> IntegralLookup<StoredMatrix>::available_basis_types() {
  std::call_once(once_register_gint_basis_types, register_gint_basis_types);

  std::vector<std::string> res(map_basis_collection_generator.size());
  std::transform(
        map_basis_collection_generator.begin(), map_basis_collection_generator.end(),
        res.begin(),
        [](const std::pair<std::string, collection_generator_type<StoredMatrix>>& kv) {
          return kv.first;
        });
  return res;
}

template <typename StoredMatrix>
IntegralLookup<StoredMatrix>::IntegralLookup(const krims::GenMap& parameters) {
  std::call_once(once_register_gint_basis_types, register_gint_basis_types);

  // TODO The orbital type is entirely ignored at the moment.
  //      We would need different integral collections for the different orbital types.
  //      Think about a way to build this into the existing system when we need it.
  const OrbitalType otype = parameters.at<OrbitalType>(IntegralLookupKeys::orbital_type,
                                                       OrbitalType::REAL_MOLECULAR);
  (void)otype;

  const auto basis_type = parameters.at<std::string>(IntegralLookupKeys::basis_type,
                                                     "<No basis_type supplied>");
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
