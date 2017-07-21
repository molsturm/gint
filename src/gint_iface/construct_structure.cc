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

#include "construct_structure.hh"

namespace gint {
namespace iface {

Structure construct_structure(long* atom_numbers, int n_atoms_an, double* coords,
                              int n_atoms_c, int three_c) {
  size_t n_atoms = static_cast<size_t>(n_atoms_an);
  assert_throw(n_atoms_an == n_atoms_c,
               krims::ExcSizeMismatch(n_atoms, static_cast<size_t>(n_atoms_c)));
  assert_throw(three_c == 3, krims::ExcSizeMismatch(static_cast<size_t>(three_c), 3ul));
  static_assert(std::is_same<real_type, double>::value, "Real type needs to be double.");

  // Build the structure object from the plain arrays:
  Structure molecule;
  molecule.reserve(n_atoms);
  for (size_t i = 0; i < n_atoms; ++i) {
    std::array<double, 3> arr_coords{
          {coords[3 * i + 0], coords[3 * i + 1], coords[3 * i + 2]}};
    assert_throw(atom_numbers[i] > 0, krims::ExcTooLarge<long>(0, atom_numbers[i]));
    molecule.emplace_back(atom_numbers[i], std::move(arr_coords));
  }
  return molecule;
}

}  // namespace iface
}  // namespace gint
