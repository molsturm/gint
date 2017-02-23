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
#include "IntegralType.hh"

namespace gint {

/** Class to uniquely identify a particular integral class.
 *
 * The class is a specialisation of a pair consisting of a string
 * (which identifies the basis type) and of an IntegralType enum
 * object.
 */
class IntegralIdentifier : std::pair<std::string, IntegralType> {
 public:
  typedef std::pair<std::string, IntegralType> base_type;
  using base_type::pair;

  /** Return the string identifying the basis type */
  std::string basis() const { return base_type::first; }

  /** Return the IntegralType object which gives the type of the integral
   * operator (overlap, kinetic, ...) */
  IntegralType integral_type() const { return base_type::second; }

  /** Return the friendly name of this integral type
   *
   * \note Is not prefixed with the basis identifier
   * */
  std::string integral_friendly_name() const { return to_friendly_name(integral_type()); }
};
}  // namespace gint
