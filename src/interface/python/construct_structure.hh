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
#include <gint/Structure.hh>

namespace gint {
namespace interface {

/** Helper function to construct a gint::Structure from a set of flat arrays. */
gint::Structure construct_structure(long* atom_numbers, int n_atoms_an, double* coords,
                                    int n_atoms_c, int three_c);

}  // namespace interface
}  // namespace gint
