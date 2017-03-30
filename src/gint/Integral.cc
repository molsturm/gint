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

#include "Integral.hh"
#include "gint/config.hh"

namespace gint {

template <typename StoredMatrix>
void Integral<StoredMatrix>::apply(
      const linalgwrap::MultiVector<const linalgwrap::MutableMemoryVector_i<scalar_type>>&
            x_in,
      linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<scalar_type>>& y_out,
      const linalgwrap::Transposed mode, const scalar_type c_this,
      const scalar_type c_y) const {
  assert_dbg(m_core_ptr != nullptr, krims::ExcInternalError());
  assert_size(x_in.n_elem(), y_out.n_elem());
  assert_size(x_in.n_vectors(), y_out.n_vectors());

  // TODO: This will go away when the new multivector interface is implemented.
  typename core_type::multivector_type x(x_in.n_elem(), x_in.n_vectors());
  typename core_type::multivector_type y(x_in.n_elem(), x_in.n_vectors());

  for (size_t i = 0; i < x_in.n_vectors(); i++)
    for (size_t j = 0; j < x_in.n_elem(); j++) {
      x(j, i) = x_in[i][j];
      y(j, i) = y_out[i][j];
    }
  m_core_ptr->apply(x, y, mode, c_this, c_y);

  for (size_t i = 0; i < x_in.n_vectors(); i++)
    for (size_t j = 0; j < x_in.n_elem(); j++) y_out[i][j] = y(j, i);
}

template <typename StoredMatrix>
void Integral<StoredMatrix>::apply_inverse(
      const linalgwrap::MultiVector<const linalgwrap::MutableMemoryVector_i<scalar_type>>&
            x_in,
      linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<scalar_type>>& y_out,
      const linalgwrap::Transposed mode, const scalar_type c_this,
      const scalar_type c_y) const {
  assert_dbg(m_core_ptr != nullptr, krims::ExcInternalError());
  assert_size(x_in.n_elem(), y_out.n_elem());
  assert_size(x_in.n_vectors(), y_out.n_vectors());

  // TODO: This will go away when the new multivector interface is implemented.
  typename core_type::multivector_type x(x_in.n_elem(), x_in.n_vectors());
  typename core_type::multivector_type y(x_in.n_elem(), x_in.n_vectors());

  for (size_t i = 0; i < x_in.n_vectors(); i++)
    for (size_t j = 0; j < x_in.n_elem(); j++) {
      x(j, i) = x_in[i][j];
      y(j, i) = y_out[i][j];
    }
  m_core_ptr->apply_inverse(x, y, mode, c_this, c_y);

  for (size_t i = 0; i < x_in.n_vectors(); i++)
    for (size_t j = 0; j < x_in.n_elem(); j++) y_out[i][j] = y(j, i);
}

// Explicitly instantiate real and complex version:
template class Integral<real_valued::stored_matrix_type>;
template class Integral<complex_valued::stored_matrix_type>;

}  // namespace gint
