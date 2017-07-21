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

#include "Element.hh"
#include "elements.hh"
#include <algorithm>

namespace gint {

bool ignore_case_equal(const std::string& a, const std::string& b) {
  return a.size() == b.size() &&
         std::equal(std::begin(a), std::end(a), std::begin(b), [](char ca, char cb) {
           return std::tolower(ca) == std::tolower(cb);
         });
}

bool is_element_symbol(const std::string& symbol) {
  return std::any_of(elements().begin(), elements().end(), [&symbol](const Element& e) {
    return ignore_case_equal(symbol, e.symbol);
  });
}

const Element& Element::by_symbol(const std::string& symbol) {
  auto res = std::find_if(
        std::begin(elements()), std::end(elements()),
        [&symbol](const Element& e) { return ignore_case_equal(symbol, e.symbol); });

  assert_throw(res != std::end(elements()),
               ExcUnknownElement("Element symbol \"" + symbol + "\" not known."));
  return *res;
}

bool is_atomic_number(unsigned int atomic_number) {
  return 0 < atomic_number && atomic_number <= elements().size();
}

const Element& Element::by_atomic_number(unsigned int atomic_number) {
  assert_throw(is_atomic_number(atomic_number),
               ExcUnknownElement("Only know atomic numbers in range [1," +
                                 std::to_string(elements().size()) + "]."));
  const auto& e = elements()[atomic_number - 1];
  assert_internal(e.atomic_number == atomic_number);
  return e;
}

}  // namespace gint
