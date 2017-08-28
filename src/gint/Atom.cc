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

Atom::Atom(const std::string& symbol, std::array<real_type, 3> coords_)
      : Atom(0, std::move(coords_)) {
  const auto& e  = Element::by_symbol(symbol);
  nuclear_charge = atomic_number = e.atomic_number;
}

std::ostream& operator<<(std::ostream& o, const Atom& atom) {
  if (is_atomic_number(atom.atomic_number)) {
    // We know this atomic number => Print element symbol
    o << Element::by_atomic_number(atom.atomic_number).symbol;
  } else {
    o << atom.atomic_number;
  }

  // If significantly different charge, print too
  if (atom.has_deviating_charge()) {
    o << " (Z=" << atom.nuclear_charge << ")";
  }

  o << "  " << atom.coords[0] << "  " << atom.coords[1] << "  " << atom.coords[2];
  return o;
}

}  // namespace gint
