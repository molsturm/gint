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
#include "gint/Structure.hh"
#include <krims/Subscribable.hh>
#include <vector>

class Shell;

namespace gint {
namespace gaussian {

/** Structure for a Gaussian basis.
 * Pretty much a slightly amended std::vector<Shell> */
struct Basis : public std::vector<Shell>, public krims::Subscribable {
  using std::vector<Shell>::vector;

  /** Construct using a molecular structure and the definition of a basis set */
  Basis(const Structure& s, const BasisSet& set);
};

}  // namespace gaussian
}  // namespace gint
