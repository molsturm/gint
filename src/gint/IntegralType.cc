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

#include "IntegralType.hh"
#include <krims/ExceptionSystem.hh>

namespace gint {
const std::string IntegralTypeKeys::nuclear_attraction = "nuclear_attraction";
const std::string IntegralTypeKeys::overlap = "overlap";
const std::string IntegralTypeKeys::kinetic = "kinetic";
const std::string IntegralTypeKeys::coulomb = "coulomb";
const std::string IntegralTypeKeys::exchange = "exchange";

#define STRING_TO_TYPE_CASE(TYPE)        \
  {                                      \
    if (key == IntegralTypeKeys::TYPE) { \
      return IntegralType::TYPE;         \
    }                                    \
  }
IntegralType to_integral_type(const std::string& key) {
  STRING_TO_TYPE_CASE(nuclear_attraction);
  STRING_TO_TYPE_CASE(overlap);
  STRING_TO_TYPE_CASE(kinetic);
  STRING_TO_TYPE_CASE(coulomb);
  STRING_TO_TYPE_CASE(exchange);

  // Else the string is not recognised
  assert_throw(false, ExcInvalidIntegralTypeKey(key));
  return IntegralType::overlap;
}

std::string to_key_string(IntegralType id) {
  switch (id) {
    case IntegralType::nuclear_attraction:
      return IntegralTypeKeys::nuclear_attraction;
    case IntegralType::overlap:
      return IntegralTypeKeys::overlap;
    case IntegralType::kinetic:
      return IntegralTypeKeys::kinetic;
    case IntegralType::coulomb:
      return IntegralTypeKeys::coulomb;
    case IntegralType::exchange:
      return IntegralTypeKeys::coulomb;
    default:
      return "Invalid ID";
  }
  return "<invalid id>";
}

std::string to_friendly_name(IntegralType id) {
  switch (id) {
    case IntegralType::nuclear_attraction:
      return "nuclear attraction";
    case IntegralType::overlap:
      return "overlap";
    case IntegralType::kinetic:
      return "kinetic";
    case IntegralType::coulomb:
      return "coulomb";
    case IntegralType::exchange:
      return "exchange";
    default:
      return "Invalid ID";
  }
  return "<invalid id>";
}

std::string to_friendly_name(const std::string& key) {
  return to_friendly_name(to_integral_type(key));
}

}  // namespace gint
