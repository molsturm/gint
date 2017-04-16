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
#include <gint/config.hh>
#include <vector>

namespace gint {
namespace gaussian {

/** Structure to represent a Gaussian shell,
 *  i.e. a set of basis functions with the same value for
 *  the azimuthal quantum number l.
 */
struct Shell {
  /** Azimuthal quantum number of the shell */
  int l;

  /** Should we use pure or Cartesian functions for the
   * spherical harmonic/angular part of the shell */
  bool pure;

  /** Contraction coefficients */
  std::vector<real_type> coefficients;

  /** Contraction exponents */
  std::vector<real_type> exponents;

  /** Origin of the shell basis functions */
  std::array<real_type, 3> origin;
};

}  // namespace gaussian
}  // namespace gint
