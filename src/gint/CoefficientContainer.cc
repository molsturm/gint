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

#include "CoefficientContainer.hh"
#include "gint/IntegralUpdateKeys.hh"

namespace gint {

template <typename StoredMatrix>
void CoefficientContainer<StoredMatrix>::update(const krims::GenMap& map) {
  const std::string& key = IntegralUpdateKeys::coefficients_occupied;
  if (!map.exists(key)) return;

  // Get coefficients as a shared pointer (having ownership)
  coeff_bo_ptr = static_cast<coefficients_ptr_type>(map.at_ptr<coefficients_type>(key));
}

// Explicitly instantiate real and complex version:
template class CoefficientContainer<real_valued::stored_matrix_type>;

// TODO complex version not yet implemented.
// template class CoefficientContainer<complex_valued::stored_matrix_type>;

}  // namespace gint
