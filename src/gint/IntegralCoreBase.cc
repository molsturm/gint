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

#include "IntegralCoreBase.hh"

namespace gint {

template <typename StoredMatrix>
void IntegralCoreBase<StoredMatrix>::apply(const const_multivector_type& x,
                                           multivector_type& y,
                                           const linalgwrap::Transposed mode,
                                           const scalar_type c_A,
                                           const scalar_type c_y) const {
  const auto& A(*this);

  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);

  // All modes are same case since we assume symmetric and real, so no
  // switching over mode.

  for (size_t i = 0; i < y.n_rows(); i++) {
    for (size_t j = 0; j < y.n_cols(); j++) {
      y(i, j) = (c_y != 0. ? c_y * y(i, j) : 0);
    }
  }

  for (size_t a = 0; a < A.n_rows(); ++a)
    for (size_t b = 0; b < A.n_cols(); ++b) {
      for (size_t q = 0; q < x.n_cols(); q++) {
        // Only implemented for symmetric matrices:
        assert_dbg(
              mode == linalgwrap::Transposed::None || std::abs(A(a, b) - A(b, a)) < 1e-14,
              krims::ExcNotImplemented());

        y(a, q) += c_A * A(a, b) * x(b, q);
      }
    }
}

template <typename StoredMatrix>
void IntegralCoreBase<StoredMatrix>::extract_block(stored_matrix_type& M,
                                                   const size_t start_row,
                                                   const size_t start_col,
                                                   const linalgwrap::Transposed mode,
                                                   const scalar_type c_A,
                                                   const scalar_type c_M) const {
  using namespace linalgwrap;

  const auto& A(*this);

  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
  assert_finite(c_A);
  assert_finite(c_M);
  // check that we do not overshoot the indices
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_greater_equal(start_row + M.n_rows(), A.n_cols());
    assert_greater_equal(start_col + M.n_cols(), A.n_rows());
  } else {
    assert_greater_equal(start_row + M.n_rows(), A.n_rows());
    assert_greater_equal(start_col + M.n_cols(), A.n_cols());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  if (c_M == 0.) {
    M.set_zero();
  } else {
    M *= c_M;
  }

  for (size_t i = start_row, i0 = 0; i < start_row + M.n_rows(); i++, i0++) {
    for (size_t j = start_col, j0 = 0; j < start_col + M.n_cols(); j++, j0++) {
      // Only implemented for symmetric matrices:
      assert_dbg(
            mode == linalgwrap::Transposed::None || std::abs(A(i, j) - A(j, i)) < 1e-14,
            krims::ExcNotImplemented());

      M(i0, j0) += c_A * A(i, j);
    }
  }
}

// Explicitly instantiate real and complex version:
template class IntegralCoreBase<real_valued::stored_matrix_type>;
template class IntegralCoreBase<complex_valued::stored_matrix_type>;

}  // namespace gint
