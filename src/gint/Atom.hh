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

namespace gint {

/** Very simple structure to describe an atom in 3D
 * A (possibly non-integer) nuclear charge and 3 coordinates. */
struct Atom {
  /** Nuclear charge of the atom
   * \note This may be set to a different value than atomic_number
   *       in order to artificially screen the charge of an atom.
   *       The mass and or the atomic number (and hence the basis
   *       functions assigned) are not changed, but the nuclear
   *       electrostatic potential is.
   */
  double nuclear_charge;

  //! Coordinates of the atom position {x,y,z}
  std::array<real_type, 3> coords;

  //! The atomic number of the atom
  unsigned int atomic_number;

  /** Construct an atom from an atomic number and a set of coordinates to place it */
  Atom(unsigned int atomic_number_, std::array<real_type, 3> coords_)
        : nuclear_charge(atomic_number_),
          coords(std::move(coords_)),
          atomic_number(atomic_number_) {}

  /** Construct an atom from an atomic symbol string and a set of coordinates */
  Atom(const std::string& symbol, std::array<real_type, 3> coords_);
};

std::ostream& operator<<(std::ostream& o, const Atom& atom);
}  // namespace gint
