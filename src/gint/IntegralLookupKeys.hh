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
#include <string>

namespace gint {

/** The struct of keys all integral libraries understand */
struct IntegralLookupKeys {
  /** The type of basis to use (Type: std::string) */
  static const std::string basis_type;

  /** The type of orbitals to use (Type: gint::OrbitalType) */
  static const std::string orbital_type;

  /** The molecular structure to model (Type: gint::Structure) */
  static const std::string structure;
};

}  // namespace gint
