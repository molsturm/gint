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
#include "Element.hh"

namespace gint {

bool ignore_case_equal(const std::string& a, const std::string& b) {
  return a.size() == b.size() &&
         std::equal(std::begin(a), std::end(a), std::begin(b), [](char ca, char cb) {
           return std::tolower(ca) == std::tolower(cb);
         });
}

Atom::Atom(const std::string& symbol, std::array<real_type, 3> coords_)
      : Atom(0.0, std::move(coords_)) {
  auto has_symbol = [&symbol](const Element& e) {
    return ignore_case_equal(symbol, e.symbol);
  };
  auto res = std::find_if(std::begin(elements()), std::end(elements()), has_symbol);

  assert_throw(res != std::end(elements()), ExcUnknownElementSymbol(symbol));
  nuclear_charge = res->atomic_number;
}

std::ostream& operator<<(std::ostream& o, const Atom& atom) {
  const int integer_charge = static_cast<int>(atom.nuclear_charge);
  if (atom.nuclear_charge - integer_charge < 1e-14 && integer_charge > 0 &&
      static_cast<size_t>(integer_charge) < elements().size()) {
    // Find element symbol and print it.

    auto has_charge = [&integer_charge](const Element& e) {
      return static_cast<unsigned short>(integer_charge) == e.atomic_number;
    };
    auto it = std::find_if(std::begin(elements()), std::end(elements()), has_charge);
    assert_dbg(it != std::end(elements()), krims::ExcInternalError());

    o << it->symbol;
  } else {
    // Print floating point charge:
    o << atom.nuclear_charge;
  }

  o << "  " << atom.coords[0] << "  " << atom.coords[1] << "  " << atom.coords[2];
  return o;
}

}  // namespace gint
