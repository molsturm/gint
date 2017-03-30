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
#include <krims/ExceptionSystem.hh>

namespace gint {

DefException1(ExcUnknownElementSymbol, std::string, << "An element with symbol \"" << arg1
                                                    << "\" is not known to gint.");

/** Very simple structure to describe an atom in 3D
 * A (possibly non-integer) nuclear charge and 3 coordinates. */
struct Atom {
  //! Nuclear charge of the atom
  double nuclear_charge;

  //! Coordinates of the atom position {x,y,z}
  std::array<real_type, 3> coords;

  Atom(double charge_, std::array<real_type, 3> coords_)
        : nuclear_charge(charge_), coords(std::move(coords_)) {}

  Atom(const std::string& symbol, std::array<real_type, 3> coords_);
};

std::ostream& operator<<(std::ostream& o, const Atom& atom);
}  // namespace gint
