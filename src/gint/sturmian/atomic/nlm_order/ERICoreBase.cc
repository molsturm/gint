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

#include "ERICoreBase.hh"
#include "gint/IntegralUpdateKeys.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

void ERICoreBase::update(const krims::GenMap& map) {
  const std::string& key = IntegralUpdateKeys::coefficients_occupied;
  if (!map.exists(key)) return;

  // Get coefficients as a shared pointer (having ownership)
  coefficients_occupied_ptr =
        static_cast<coefficients_ptr_type>(map.at_ptr<coefficients_type>(key));

  // We will contract the coefficient row index over the number of
  // basis functions.
  if (coefficients_occupied_ptr->n_vectors() == 0) return;
  assert_size(coefficients_occupied_ptr->n_elem(), system().n_bas());
}

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
