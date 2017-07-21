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
#include "BasisSet.hh"
#include "Shell.hh"
#include <krims/Subscribable.hh>
#include <vector>

namespace gint {

#ifndef SWIG
// Forward-declare Structure class
class Structure;
#endif  // SWIG

namespace gaussian {

/** Structure for a Gaussian basis.
 * Pretty much a slightly amended std::vector<Shell> */
struct Basis : public std::vector<Shell>, public krims::Subscribable {
#ifndef SWIG
  using std::vector<Shell>::vector;

  /** Construct using a molecular structure and the definition of a basis set */
  Basis(const Structure& s, const BasisSet& set);
#endif  // SWIG

  /** Return the number of shells in the basis */
  size_t n_shells() const { return size(); }

  /** Default constructor */
  Basis() {}
};

#if SWIG
/* clang-format off */
%extend Basis {
  const gint::gaussian::Shell& shell(size_t i) const { return $self->operator[](i); }
}
/* clang-format on */
#endif  // SWIG

}  // namespace gaussian
}  // namespace gint
