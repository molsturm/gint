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

#include "construct_gaussian_basis.hh"
#include <gint/gaussian/BasisSet.hh>

namespace gint {
namespace iface {

gint::gaussian::Basis construct_gaussian_basis(long* atom_numbers, int n_atoms_an,
                                               double* coords, int n_atoms_c, int three_c,
                                               const std::string& basis_name) {
  gint::Structure molecule =
        construct_structure(atom_numbers, n_atoms_an, coords, n_atoms_c, three_c);
  return gint::gaussian::Basis(molecule, gaussian::lookup_basisset(basis_name));
}

}  // namespace iface
}  // namespace gint
