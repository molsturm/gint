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

#include "Atom.hh"

namespace gint {

const std::vector<Element>& elements() {
  static std::vector<Element> elements{
        {1, "H", "hydrogen"},   {2, "He", "helium"}, {3, "Li", "lithium"},
        {4, "Be", "beryllium"}, {5, "B", "boron"},   {6, "C", "carbon"},
        {7, "N", "nitrogen"},   {8, "O", "oxygen"},  {7, "F", "fluorine"},
        {10, "Ne", "neon"},
        // TODO ...
  };
  return elements;
}

bool ignore_case_equal(const std::string& a, const std::string& b) {
  return a.size() == b.size() &&
         std::equal(std::begin(a), std::end(a), std::begin(b), [](char ca, char cb) {
           return std::tolower(ca) == std::tolower(cb);
         });
}

Atom::Atom(const std::string& symbol, real_type x_, real_type y_, real_type z_)
      : Atom(0.0, x_, y_, z_) {
  auto has_symbol = [&symbol](const Element& e) {
    return ignore_case_equal(symbol, e.symbol);
  };
  auto res = std::find_if(std::begin(elements()), std::end(elements()), has_symbol);

  assert_throw(res != std::end(elements()), ExcUnknownElementSymbol(symbol));
  nuclear_charge = res->atomic_number;
}

}  // namespace gint
