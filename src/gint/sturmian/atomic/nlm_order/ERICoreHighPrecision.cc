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

#include "ERICoreHighPrecision.hh"
#ifdef GINT_HAVE_STURMINT

// For explicit instatiations
#include <sturmint/atomic/cs_naive/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

template <typename RepulsionCalculator>
void ERICoreHighPrecision<RepulsionCalculator>::apply(const const_multivector_type& x,
                                                      multivector_type& y,
                                                      const lazyten::Transposed mode,
                                                      const scalar_type c_A,
                                                      const scalar_type c_y) const {
  // TODO Conceptionally this is code duplication with the apply function
  // in gint/IntegralCoreBase.cc, but just performed at elevated precision.
  //
  // If the elevated precision is not need, one might be able to remove this,
  // else one could merge it into the above place.

  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), base_type::n_cols());
  assert_size(y.n_rows(), base_type::n_rows());
  assert_sufficiently_tested(mode != lazyten::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  const size_t norb = x.n_rows(), n_vectors = x.n_cols();

  for (size_t i = 0; i < y.n_rows(); i++) {
    for (size_t j = 0; j < y.n_cols(); j++) {
      y(i, j) = (c_y != 0 ? c_y * y(i, j) : 0);
    }
  }

  // Work internally in highest available precision
  std::vector<working_scalar_type> ysum(n_vectors);
  for (size_t a = 0; a < norb; a++) {
    std::fill(std::begin(ysum), std::end(ysum), 0);
    for (size_t b = 0; b < norb; b++) {
      for (size_t q = 0; q < n_vectors; q++) {
        ysum[q] += c_A * base_type::compute_jk_element(a, b) * x(b, q);
      }
    }
    for (size_t q = 0; q < n_vectors; q++) y(a, q) += ysum[q];
  }
}

#define INSTANTIATE(CLASS)                                                     \
  template void ERICoreHighPrecision<CLASS>::apply(                            \
        const const_multivector_type&, multivector_type&, lazyten::Transposed, \
        const scalar_type, const scalar_type) const

INSTANTIATE(sturmint::atomic::cs_naive::Atomic);

#undef INSTANTIATE

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
