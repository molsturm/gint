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
#include <krims/ExceptionSystem.hh>
#include <map>
#include <vector>

namespace gint {
namespace gaussian {
struct Shell;

enum class DefaultAngularFunctions {
  /** Use pure spherical harmonics for representing the angular momentum */
  Pure = 0,

  /** Use cartesian functions for representing the angular momentum */
  Cartesian = 1,

  /** Use cartesian functions only for d, but pure spherical harmonics for higher
     angular momentum */
  CartesianForD = 3,
};

DefException3(ExcNoBasisForAtom, unsigned int, std::string, std::string,
              << "Could not find atom number " << arg1 << " in basis set \"" << arg2
              << "\", which was read from file \"" << arg3 << "\".");

/** A structure which should represent a parsed Gaussian basis set
 *  (e.g. something read in from a file */
struct BasisSet {
  /** Name of the basis set (e.g. "sto-3g", "cc-pvtz") */
  std::string name{"<unknown>"};

  /** Path to the file from which this basis set was read. */
  std::string filename{"<unknown>"};

  /** Obtain the shells which should be used to model a particular atom
   *  in this basis set.
   *
   * \throws ExcNoBasisForAtom if this basis does not contain any shells
   * for a particular atom.
   */
  const std::vector<Shell>& shells_for_atom(unsigned int atomic_number) const {
    auto it = atomic_number_to_shells.find(atomic_number);
    assert_throw(it != std::end(atomic_number_to_shells),
                 ExcNoBasisForAtom(atomic_number, name, filename));
    return it->second;
  }

  /** The mapping from the atomic number Z to the list of
   *  Shells which should be used for modelling it.
   */
  std::map<unsigned int, std::vector<Shell>> atomic_number_to_shells;
};

/** Lookup basis set by name.
 *
 * Tries to find the basis with the data distributed with gint.
 * The input string is normalised before use.
 *
 * \note This function honours the convention imposed by the function
 *       lookup_default_angular_functions when it comes to the kind of
 *       angular momentum functions used (Cartesian or Pure).
 */
BasisSet lookup_basisset(const std::string& name);

/** Determine basis set options from the name
 *
 * Some basis sets (like 6-31G*) enforce the use of Cartesian functions for the
 * d orbitals. This sets the relevant options in those cases
 * See http://gaussian.com/basissets
 */
DefaultAngularFunctions lookup_default_angular_functions(const std::string& name);

}  // namespace gaussian
}  // namespace gint
