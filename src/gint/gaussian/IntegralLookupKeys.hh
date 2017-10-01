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
#include "gint/IntegralLookupKeys.hh"
#include <string>

namespace gint {
namespace gaussian {

/** The struct of keys all gaussian integral library integral collections understand */
struct IntegralLookupKeys : public gint::IntegralLookupKeys {
  /** The basis set to use (Type: std::string)
   *  Either the name of a basis set, which is looked up on the filesystem,
   *  or a (realative) path to a basis set to use.
   *
   *  The actual basis, which is used is constructed from the obtained basis set
   *  using the structure to model.
   */
  static const std::string basis_set_name;

  /** The basis to use (Type: gint::gaussian::Basis)
   *  Takes preference over the basis_set supplied */
  static const std::string basis;
};

}  // namespace gaussian
}  // namespace gint
