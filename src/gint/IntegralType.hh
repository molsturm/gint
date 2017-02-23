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

#pragma once
#include <krims/ExceptionSystem.hh>
#include <limits>
#include <string>

namespace gint {

/** Enum describing the types of integrals available */
enum class IntegralType {
  /** Nuclear attraction integral */
  nuclear_attraction,

  /** Overlap integral */
  overlap,

  /** Kinetic integral */
  kinetic,

  /** Coulomb integral */
  coulomb,

  /** Exchange integral */
  exchange,
};

/** Struct which holds the static id strings used to refer to the individual
 * integral types */
struct IntegralTypeKeys {
  /** Name string used for nuclear attraction integrals */
  const static std::string nuclear_attraction;

  /** Name string used for overlap integrals */
  const static std::string overlap;

  /** Name string used for kinetic integrals */
  const static std::string kinetic;

  /** Name string used for coulomb integrals */
  const static std::string coulomb;

  /** Name string used for exchange integrals */
  const static std::string exchange;
};

/** Return the key string of the IntegralTypeKeys struct which matches the IntegralType */
std::string to_key_string(IntegralType id);

//@{
/** Return the friendly name of the integral type */
std::string to_friendly_name(IntegralType id);
std::string to_friendly_name(const std::string& key);
//@}

DefException1(ExcInvalidIntegralTypeKey, std::string, << "The integral type key string \""
                                                      << arg1 << "\" is unknown.");

/** Return the IntegralType which corresponds to the key string
 *
 * \throws ExcInvalidIntegralTypeKey in case the string is not a valid key string.
 */
IntegralType to_integral_type(const std::string& key);

}  // namespace gint
