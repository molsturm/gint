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

#include "ERITensor.hh"
#ifdef GINT_HAVE_STURMINT

// For explicit instatiations
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <sturmint/atomic/cs_naive/cs_atomic.hh>
#include <sturmint/atomic/cs_reference/cs_atomic.hh>
#include <sturmint/atomic/cs_reference_pc/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

template <typename RepulsionCalculatior>
void ERITensor<RepulsionCalculatior>::compute_kernel(
      const std::array<krims::Range<size_t>, 4>& block, kernel_type kernel) const {

  // Do the computation in 8 x 8 x 8 x 8 block
  // TODO This value is arbitrarily chosen!
  const size_t bs = 8;

  // Get a buffer on the stack
  std::array<scalar_type, bs * bs * bs * bs> buffer;

  for (size_t ab = block[0].lower_bound(); ab < block[0].upper_bound(); ab += bs) {
    for (size_t bb = block[1].lower_bound(); bb < block[1].upper_bound(); bb += bs) {
      for (size_t cb = block[2].lower_bound(); cb < block[2].upper_bound(); cb += bs) {
        for (size_t db = block[3].lower_bound(); db < block[3].upper_bound(); db += bs) {
          const krims::Range<size_t> a_range{
                {ab, std::min(block[0].upper_bound(), ab + bs)}};
          const krims::Range<size_t> b_range{
                {bb, std::min(block[1].upper_bound(), bb + bs)}};
          const krims::Range<size_t> c_range{
                {cb, std::min(block[2].upper_bound(), cb + bs)}};
          const krims::Range<size_t> d_range{
                {db, std::min(block[3].upper_bound(), db + bs)}};

          // Fill buffer
          auto it = buffer.begin();
          for (auto a : a_range) {
            for (auto b : b_range) {
              for (auto c : c_range) {
                for (size_t d = d_range.lower_bound(); d < d_range.upper_bound();
                     ++d, ++it) {
                  std::array<int, 4> idcs{{static_cast<int>(a), static_cast<int>(b),
                                           static_cast<int>(c), static_cast<int>(d)}};
                  *it = m_system_ptr->k *
                        static_cast<scalar_type>(m_calculator_ptr->repulsion(idcs));
                }  // d
              }    // c
            }      // b
          }        // a

          kernel({{a_range, b_range, c_range, d_range}}, buffer.data());

        }  // db
      }    // cb
    }      // bb
  }        // ab
}

#define INSTANTIATE(CLASS)                        \
  template void ERITensor<CLASS>::compute_kernel( \
        const std::array<krims::Range<size_t>, 4>&, kernel_type) const

INSTANTIATE(sturmint::atomic::cs_dummy::Atomic);
INSTANTIATE(sturmint::atomic::cs_naive::Atomic);
INSTANTIATE(sturmint::atomic::cs_reference::Atomic);
INSTANTIATE(sturmint::atomic::cs_reference_pc::Atomic);

#undef INSTANTIATE

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
