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
#include "Atom.hh"
#include <krims/Subscribable.hh>

namespace gint {

/** Data structure describing the geometry of a molecular structure as
 *  a std::vector of atoms */
class Structure : public std::vector<Atom>, public krims::Subscribable {
 public:
  typedef std::vector<Atom> base_type;

  // Use constructors from std::vector
  using std::vector<Atom>::vector;

  /** Return the number of atoms in this molecular structure */
  size_t n_atoms() const { return base_type::size(); }

  /** Return the total charge of all atoms in this molecular structure */
  double total_charge() const {
    double ret{0};
    for (const auto& a : (*this)) ret += a.nuclear_charge;
    return ret;
  }
};

std::ostream& operator<<(std::ostream& o, const Structure& st);

}  // namespace gint
