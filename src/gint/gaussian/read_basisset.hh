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
#include <istream>
#include <krims/ExceptionSystem.hh>

namespace gint {
namespace gaussian {

struct BasisSet;

/** Enum allowing to multiplex between the different basis set file formats gint
 * understands. */
enum class BasisSetFileFormat {
  /** Gaussian 94 basis format, see http://gaussian.com/gen/ for details */
  Gaussian94
};

DefException1(ExcInvalidBasisSetFile, std::string, << arg1);

/** Read a basis set from a file and return a BasisSet structure
 *
 * \throws ExcInvalidBasisSetFile if the file cannot be parsed.
 * */
BasisSet read_basisset(std::istream& in, BasisSetFileFormat fmt);

}  // namespace gaussian
}  // namespace gint
