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

#include "ERITensor_i.hh"
#include "gint/config.hh"

namespace gint {

template <typename Scalar>
void ERITensor_i<Scalar>::contract_with(const iface_multivector_type& c_wa,
                                        const iface_multivector_type& c_xb,
                                        const iface_multivector_type& c_yc,
                                        const iface_multivector_type& c_zd,
                                        std::vector<Scalar>& out) const {
  assert_size(c_wa.n_vectors() * c_xb.n_vectors() * c_yc.n_vectors() * c_zd.n_vectors(),
              out.size());

  if (out.size() == 0) return;  // Nothing to be done

  assert_size(c_wa.n_elem(), n_bas());
  assert_size(c_xb.n_elem(), n_bas());
  assert_size(c_yc.n_elem(), n_bas());
  assert_size(c_zd.n_elem(), n_bas());

  auto contract_kernel = [&out, &c_wa, &c_xb, &c_yc, &c_zd](
        const std::array<krims::Range<size_t>, 4>& batch, const Scalar* values) {
    const size_t n_w = c_wa.n_vectors();
    const size_t n_x = c_xb.n_vectors();
    const size_t n_y = c_yc.n_vectors();
    const size_t n_z = c_zd.n_vectors();

#ifdef DEBUG
    const size_t values_size =
          batch[0].length() * batch[1].length() * batch[2].length() * batch[3].length();
#endif
    for (size_t w = 0, i_wxyz = 0; w < n_w; ++w) {
      for (size_t x = 0; x < n_x; ++x) {
        for (size_t y = 0; y < n_y; ++y) {
          for (size_t z = 0; z < n_z; ++z, ++i_wxyz) {

            const Scalar* it = values;
            for (auto a : batch[0]) {
              for (auto b : batch[1]) {
                for (auto c : batch[2]) {
                  for (auto d = batch[3].lower_bound(); d < batch[3].upper_bound();
                       ++d, ++it) {
                    assert_internal(static_cast<size_t>(it - values) < values_size);
                    assert_internal(i_wxyz < out.size());
                    out[i_wxyz] +=
                          c_wa[w][a] * c_xb[x][b] * c_yc[y][c] * c_zd[z][d] * *it;
                  }  // d
                }    // c
              }      // b
            }        // a

          }  // z
        }    // y
      }      // x
    }        // w
  };
  compute_kernel(std::move(contract_kernel));
}

template <typename Scalar>
void ERITensor_i<Scalar>::extract_block(const std::array<krims::Range<size_t>, 4>& block,
                                        std::vector<Scalar>& out) const {
  assert_size(out.size(), block[0].length() * block[1].length() * block[2].length() *
                                block[3].length());
  assert_greater(block[0].back(), n_bas());
  assert_greater(block[1].back(), n_bas());
  assert_greater(block[2].back(), n_bas());
  assert_greater(block[3].back(), n_bas());

  auto extract_kernel = [&out, &block](const std::array<krims::Range<size_t>, 4>& batch,
                                       const Scalar* values) {
#ifdef DEBUG
    const size_t values_size =
          batch[0].length() * batch[1].length() * batch[2].length() * batch[3].length();
#endif

    const Scalar* it = values;
    for (auto a : batch[0]) {
      for (auto b : batch[1]) {
        for (auto c : batch[2]) {
          const size_t abc = ((a * block[0].length() + b) * block[1].length() + c);

          for (auto d = batch[3].lower_bound(); d < batch[3].upper_bound(); ++d, ++it) {
            const size_t abcd = abc * block[2].length() + d;
            assert_internal(static_cast<size_t>(it - values) < values_size);
            assert_internal(abcd < out.size());
            out[abcd] = *it;
          }  // d
        }    // c
      }      // b
    }        // a
  };

  compute_kernel(block, std::move(extract_kernel));
}

// Explicit instatiations
template class ERITensor_i<real_type>;
template class ERITensor_i<complex_type>;

}  // namespace gint
