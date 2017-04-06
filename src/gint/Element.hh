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
#include <array>
#include <krims/ExceptionSystem.hh>
#include <string>

namespace gint {

DefException1(ExcUnknownElement, std::string, << "Element unknown to gint: " << arg1);

/** Structure for all elements of the periodic table */
struct Element {
  unsigned int atomic_number;
  std::string symbol;
  std::string name;

  /** Look up an element by element symbol.
   * \throws ExcUnknownElement if the symbol is not known */
  static const Element& by_symbol(const std::string& symbol);

  /** Look up an element by atomic number.
   * \throws ExcUnknownElement if the atomic number is not known */
  static const Element& by_atomic_number(unsigned int atomic_number);
};

/** Get the list of all elements known to gint
 * No particular order should be assumed here.
 * \note Use Element::by_symbol or Element::by_atomic_number
 * to lookup an element.
 */
const std::array<Element, 118>& elements();

/** Check whether the provided number is a valid atomic number */
bool is_atomic_number(unsigned int atomic_number);

/** Check whether the provided string is a valid element symbol */
bool is_element_symbol(const std::string& symbol);

}  // namespace gint
