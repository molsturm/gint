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

#ifdef SWIG
/* clang-format off */
%apply(int* DIM1, double** ARGOUTVIEW_ARRAY1){(int* n_coeff, double** coeff)};
%apply(int* DIM1, double** ARGOUTVIEW_ARRAY1){(int* n_exp,   double** exp  )};
%apply(int* DIM1, double** ARGOUTVIEW_ARRAY1){(int* n_coord, double** coord)};
/* clang-format on */
#endif  // SWIG

namespace gint {
namespace gaussian {

#ifndef SWIG
static_assert(std::is_same<double, real_type>::value,
              "Currently real_type needs to be double for SWIG interoperability.");
#endif

/** Structure to represent a Gaussian shell,
 *  i.e. a set of basis functions with the same value for
 *  the azimuthal quantum number l.
 */
struct Shell {
  /** Azimuthal quantum number of the shell */
  int l;

  /** Should we use pure or Cartesian functions for the
   * spherical harmonic/angular part of the shell
   *
   * \note
   * For angular momentum quantum numbers less than
   * 2 (i.e. s and p shells) the value of pure should not
   * make a difference. By convention it should be false
   * in these cases.
   **/
  bool pure;

#ifndef SWIG
  //@{
  /** Return the number of basis functions in this shell */
  int size() const { return pure ? (2 * l + 1) : (l + 1) * (l + 2) / 2; }
  int n_bas() const { return size(); }
  //@}

  //@{
  /** Return the number of primitives, i.e. the number of basis functions which are
   * contracted together */
  size_t n_primitives() const { return coefficients.size(); }
  size_t n_contracted() const { return n_primitives(); }
  //@}

  /** Contraction coefficients */
  std::vector<real_type> coefficients;

  /** Contraction exponents */
  std::vector<real_type> exponents;

  /** Origin of the shell basis functions */
  std::array<real_type, 3> origin;
#endif  // SWIG
};

#if SWIG
/* clang-format off */
%extend Shell {
  void coefficients(int* n_coeff, double** coeff) {
    *n_coeff = $self->coefficients.size();
    *coeff   = $self->coefficients.data();
  }

  void exponents(int* n_exp, double** exp) {
    *n_exp = $self->exponents.size();
    *exp   = $self->exponents.data();
  }

  void origin(int* n_coord, double** coord) {
    *n_coord = 3;
    *coord   = $self->origin.data();
  }
}
/* clang-format on */
#endif  // SWIG

}  // namespace gaussian
}  // namespace gint
