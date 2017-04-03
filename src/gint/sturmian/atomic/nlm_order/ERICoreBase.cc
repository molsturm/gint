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

// For explicit instatiation below ...
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <sturmint/atomic/cs_reference/cs_atomic.hh>
#include <sturmint/atomic/cs_reference_pc/cs_atomic.hh>

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
  assert_size(coefficients_occupied_ptr->n_elem(), n_bas());
}

template <typename Calculator>
scalar_type ERICoreBase::compute_jk_element(const Calculator& calculator, size_t a,
                                            size_t b) const {
  assert_greater(a, n_rows());
  assert_greater(b, n_cols());
  assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

  // Compute density matrix:
  const auto density = linalgwrap::outer_prod_sum(coeff_bo(), coeff_bo());

  scalar_type sum = 0;
  for (size_t c = 0; c < n_bas(); c++) {
    // Swap a and c if computing exchange
    const bool exchange = type() == IntegralType::exchange;
    const size_t A = exchange ? c : a;
    const size_t C = exchange ? a : c;

    for (size_t d = 0; d < n_bas(); d++)
      sum += calculator.repulsion(A, b, C, d) * density(c, d);
  }
  return system().k * sum;
}

#define INSTANTIATE_COMPUTE_JK_ELEMENT(CLASS)                                          \
  template scalar_type ERICoreBase::compute_jk_element(const sturmint::atomic::CLASS&, \
                                                       size_t, size_t) const

INSTANTIATE_COMPUTE_JK_ELEMENT(cs_dummy::Atomic);
INSTANTIATE_COMPUTE_JK_ELEMENT(cs_reference::Atomic);
INSTANTIATE_COMPUTE_JK_ELEMENT(cs_reference_pc::Atomic);

#undef INSTANTIATE_COMPUTE_JK_ELEMENT

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
