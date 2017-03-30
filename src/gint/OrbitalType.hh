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

namespace gint {

/** Enum class describing the type of orbital requested */
enum class OrbitalType {
  /** Real atomic orbitals */
  REAL_ATOMIC,

  /** Complex atomic orbitals */
  COMPLEX_ATOMIC,

  /** Real molecular orbitals */
  REAL_MOLECULAR,

  /** Complex molecular orbitals */
  COMPLEX_MOLECULAR
};

}  // namespace gint
