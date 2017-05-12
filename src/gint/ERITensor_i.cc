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

template <size_t N>
struct MakeAbsolute {
  size_t operator()(std::array<size_t, N> idcs) {
    auto itidx = idcs.begin();
    auto itshape = shape.begin();

    size_t accu = 0;
    for (; itidx < idcs.end(); ++itidx, ++itshape) {
      accu = accu * (*itshape) + (*itidx);
    }
    return accu;
  }
  std::array<size_t, N> shape;

  explicit MakeAbsolute(std::array<size_t, N> shape_) : shape(std::move(shape_)) {}
};

template <typename Scalar>
void ERITensor_i<Scalar>::contract_with(const iface_multivector_type& c_wa,
                                        const iface_multivector_type& c_xb,
                                        const iface_multivector_type& c_yc,
                                        const iface_multivector_type& c_zd,
                                        std::vector<Scalar>& out) const {
  // The contraction proceeds in four steps, which are each O(n_bas^4) in memory
  // and O(n_bas^4 n_orb) in computational time.
  //
  // So essentially we do
  //  \f[    I1_{wbcd} = \sum_a (ab|cd) C_{wa}  \f]
  // then
  //  \f[    I2_{wxcd} = \sum_a I1_{wbcd} C_{xb}  \f]
  // and
  //  \f[    I3_{wxyd} = \sum_a I2_{wxcd} C_{yc}  \f]
  // and finally
  //  \f[    out_{wxyz} = \sum_a I3_{wxyd} C_{zd}  \f]

  assert_size(c_wa.n_vectors() * c_xb.n_vectors() * c_yc.n_vectors() * c_zd.n_vectors(),
              out.size());

  if (out.size() == 0) return;  // Nothing to be done

  assert_size(c_wa.n_elem(), n_bas());
  assert_size(c_xb.n_elem(), n_bas());
  assert_size(c_yc.n_elem(), n_bas());
  assert_size(c_zd.n_elem(), n_bas());

  const size_t n_w = c_wa.n_vectors();
  const size_t n_x = c_xb.n_vectors();
  const size_t n_y = c_yc.n_vectors();
  const size_t n_z = c_zd.n_vectors();

  // Index transformation helpers
  MakeAbsolute<4> absidx_wbcd({{n_w, n_bas(), n_bas(), n_bas()}});
  MakeAbsolute<4> absidx_wxcd({{n_w, n_x, n_bas(), n_bas()}});
  MakeAbsolute<4> absidx_wxyd({{n_w, n_x, n_y, n_bas()}});
  MakeAbsolute<4> absidx_wxyz({{n_w, n_x, n_y, n_z}});

  //
  // Build intermediate 1
  //
  std::vector<Scalar> i1_wbcd(n_w * n_bas() * n_bas() * n_bas());

  auto contract_kernel = [&i1_wbcd, &c_wa, &absidx_wbcd](
        const std::array<krims::Range<size_t>, 4>& batch, const Scalar* values) {
#ifdef DEBUG
    const size_t values_size =
          batch[0].length() * batch[1].length() * batch[2].length() * batch[3].length();
#endif

    for (size_t w = 0; w < c_wa.n_vectors(); ++w) {
      size_t abcd = 0;
      for (auto a : batch[0]) {

        for (auto b : batch[1]) {
          for (auto c : batch[2]) {
            for (size_t d = batch[3].lower_bound(); d < batch[3].upper_bound();
                 ++d, ++abcd) {

              const size_t wbcd = absidx_wbcd({{w, b, c, d}});
              assert_internal(abcd < values_size);
              assert_internal(wbcd < i1_wbcd.size());
              i1_wbcd[wbcd] += c_wa[w][a] * values[abcd];
            }  // d
          }    // c
        }      // b

      }  // a
    }    // w
  };
  compute_kernel(std::move(contract_kernel));

  //
  // Build intermediate 2
  //
  auto& i2_wxcd(out);  // Re-use memory for out
  i2_wxcd.resize(n_w * n_x * n_bas() * n_bas());
  std::fill(i2_wxcd.begin(), i2_wxcd.end(), 0);

  for (size_t w = 0; w < n_w; ++w) {

    for (size_t x = 0; x < n_x; ++x) {
      for (size_t b = 0; b < n_bas(); ++b) {

        for (size_t c = 0; c < n_bas(); ++c) {
          for (size_t d = 0; d < n_bas(); ++d) {
            const size_t wbcd = absidx_wbcd({{w, b, c, d}});
            const size_t wxcd = absidx_wxcd({{w, x, c, d}});
            assert_internal(wbcd < i1_wbcd.size());
            assert_internal(wxcd < i2_wxcd.size());

            i2_wxcd[wxcd] += i1_wbcd[wbcd] * c_xb[x][b];
          }  // d
        }    // c

      }  // b
    }    // x

  }  // w

  //
  // Build intermediate 3
  //
  auto& i3_wxyd = i1_wbcd;  // Re-use memory for i3
  i3_wxyd.resize(n_w * n_x * n_y * n_bas());
  std::fill(i3_wxyd.begin(), i3_wxyd.end(), 0);

  for (size_t w = 0; w < n_w; ++w) {
    for (size_t x = 0; x < n_x; ++x) {

      for (size_t y = 0; y < n_y; ++y) {
        for (size_t c = 0; c < n_bas(); ++c) {

          for (size_t d = 0; d < n_bas(); ++d) {
            const size_t wxcd = absidx_wxcd({{w, x, c, d}});
            const size_t wxyd = absidx_wxyd({{w, x, y, d}});
            assert_internal(wxcd < i2_wxcd.size());
            assert_internal(wxyd < i3_wxyd.size());

            i3_wxyd[wxyd] += i2_wxcd[wxcd] * c_yc[y][c];
          }  // d

        }  // c
      }    // y

    }  // x
  }    // w

  //
  // Build intermediate 4
  //
  out.resize(n_w * n_x * n_y * n_z);
  std::fill(out.begin(), out.end(), 0);

  for (size_t w = 0; w < n_w; ++w) {
    for (size_t x = 0; x < n_x; ++x) {
      for (size_t y = 0; y < n_y; ++y) {

        for (size_t z = 0; z < n_z; ++z) {
          for (size_t d = 0; d < n_bas(); ++d) {
            const size_t wxyd = absidx_wxyd({{w, x, y, d}});
            const size_t wxyz = absidx_wxyz({{w, x, y, z}});
            assert_internal(wxyd < i3_wxyd.size());
            assert_internal(wxyz < out.size());

            out[wxyz] += i3_wxyd[wxyd] * c_zd[z][d];
          }  // d
        }    // z

      }  // y
    }    // x
  }      // w
}

template <typename Scalar>
void ERITensor_i<Scalar>::extract_block(const std::array<krims::Range<size_t>, 4>& block,
                                        std::vector<Scalar>& out) const {
  using krims::intersection;
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

    // Get the intersection of the range we have in the current batch
    // and the block of the tensor we want to extract.
    // So this is now the range of indices we care about in this current
    // execution of the lambda.
    const auto a_range = intersection(batch[0], block[0]);
    const auto b_range = intersection(batch[1], block[1]);
    const auto c_range = intersection(batch[2], block[2]);
    const auto d_range = intersection(batch[3], block[3]);

    // Get absolute indices for indexing into the current batch
    // or into the block of interest
    MakeAbsolute<4> absidx_blk(
          {{block[0].length(), block[1].length(), block[2].length(), block[3].length()}});
    MakeAbsolute<4> absidx_bch(
          {{batch[0].length(), batch[1].length(), batch[2].length(), batch[3].length()}});

    const Scalar* it = values;
    for (auto a : a_range) {
      const size_t a_bch = a - batch[0].front();
      const size_t a_blk = a - block[0].front();

      for (auto b : b_range) {
        const size_t b_bch = b - batch[1].front();
        const size_t b_blk = b - block[1].front();

        for (auto c : c_range) {
          const size_t c_bch = c - batch[2].front();
          const size_t c_blk = c - block[2].front();

          for (auto d : d_range) {
            const size_t d_bch = d - batch[3].front();
            const size_t d_blk = d - block[3].front();

            const size_t idx_blk = absidx_blk({{a_blk, b_blk, c_blk, d_blk}});
            const size_t idx_bch = absidx_bch({{a_bch, b_bch, c_bch, d_bch}});

            assert_internal(idx_bch < values_size);
            assert_internal(idx_blk < out.size());
            out[idx_blk] = it[idx_bch];
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
