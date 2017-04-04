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

#include "ERICore.hh"

// For explicit instatiations of compute_jk_element
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <sturmint/atomic/cs_naive/cs_atomic.hh>
#include <sturmint/atomic/cs_reference/cs_atomic.hh>
#include <sturmint/atomic/cs_reference_pc/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {
typedef sturmint::real_type working_scalar_type;

template <typename RepulsionCalculator>
working_scalar_type ERICore<RepulsionCalculator>::compute_jk_element(size_t a,
                                                                     size_t b) const {
  assert_greater(a, n_rows());
  assert_greater(b, n_cols());

  // TODO The density matrix expression: Maybe bump to working_scalar_type as well!
  const auto dens = compute_density_matrix();

  /*
   *  working_scalar_type density_cd{0};
      for (size_t p = 0; p < coeff_bo().n_vectors(); p++)
        density_cd += coeff_bo()[p][c] * coeff_bo()[p][d];
  */

  working_scalar_type sum = 0;
  for (size_t c = 0; c < n_bas(); c++) {
    // Swap a and c if computing exchange
    const bool exchange = type() == IntegralType::exchange;
    const size_t A = exchange ? c : a;
    const size_t C = exchange ? a : c;

    for (size_t d = 0; d < n_bas(); d++) {
      std::array<int, 4> idcs{{static_cast<int>(A), static_cast<int>(b),
                               static_cast<int>(C), static_cast<int>(d)}};
      sum += m_calculator_ptr->repulsion(idcs) * dens(c, d);
    }
  }
  return system().k * sum;
}

#define INSTANTIATE(CLASS) \
  template working_scalar_type ERICore<CLASS>::compute_jk_element(size_t, size_t) const

INSTANTIATE(sturmint::atomic::cs_dummy::Atomic);
INSTANTIATE(sturmint::atomic::cs_naive::Atomic);
INSTANTIATE(sturmint::atomic::cs_reference::Atomic);
INSTANTIATE(sturmint::atomic::cs_reference_pc::Atomic);

#undef INSTANTIATE

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
