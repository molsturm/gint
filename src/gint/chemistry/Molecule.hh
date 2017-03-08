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
#include "gint/config.hh"
#include <krims/Subscribable.hh>

namespace gint {

/** Very simple structure to describe an atom in 3D
 * A (possibly non-integer) nuclear charge and 3 coordinates.
 */
struct Atom {
  float nuclear_charge;
  real_type x, y, z;
};

/** Data structure describing a molecule as a std::vector of atoms */
class Molecule : public std::vector<Atom>, public krims::Subscribable {
 public:
  using std::vector<Atom>::vector;

  typedef std::vector<Atom> base_type;

  /** Return the number of atoms in this molecule */
  size_t n_atoms() const { return base_type::size(); }
};

}  // namespace gint
