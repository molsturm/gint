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

#pragma once
#include <linalgwrap/MultiVector.hh>

namespace gint {

template <typename Scalar>
class EriTensor_i {
 public:
  typedef const linalgwrap::MutableMemoryVector_i<Scalar> iface_vector_type;
  typedef linalgwrap::MultiVector<iface_vector_type> iface_multivector_type;

  /** Contract the electron repulsion tensor, which is implicitly represented
   *  by this structure, with a set of multivectors.
   *
   * Essentially it carries out the tensor contraction
   * \f[
   *   J_{wx,yz} = \sum_{b1,b2,b3,b4} c^T_{w,b1} c_{w,b2}, c^T_{w,b3}, c_z,b4}
   * J_{b1b2,b3b4}
   * \f]
   * where \f$ J_{b1b2,b3b4} = (b1 b2| b3 b4) \f$, i.e. the first two and the last
   * two sit on the same centre.
   *
   * The resulting rank 4 tensor is written in linearised fashion into the output
   * vector, such that later indices run faster (C-like convention).
   */
  template <typename Vector>
  void contract_with(const linalgwrap::MultiVector<Vector>& c_wb1,
                     const linalgwrap::MultiVector<Vector>& c_xb2,
                     const linalgwrap::MultiVector<Vector>& c_yb3,
                     const linalgwrap::MultiVector<Vector>& c_zb4,
                     std::vector<Scalar>& out) const {

    using namespace linalgwrap;
    iface_multivector_type c_wb1_wrap(c_wb1);
    iface_multivector_type c_xb2_wrap(c_xb2);
    iface_multivector_type c_yb3_wrap(c_yb3);
    iface_multivector_type c_zb4_wrap(c_zb4);

    contract_with(c_wb1_wrap, c_xb2_wrap, c_yb3_wrap, c_zb4_wrap, out);
  }

  /** Contract the electron repulsion tensor, which is implicitly represented
   *  by this structure, with a set of multivectors.
   *
   * Essentially it carries out the tensor contraction
   * \f[
   *   J_{wx,yz} = \sum_{b1,b2,b3,b4} c^T_{w,b1} c_{w,b2}, c^T_{w,b3}, c_z,b4}
   * J_{b1b2,b3b4}
   * \f]
   * where \f$ J_{b1b2,b3b4} = (b1 b2| b3 b4) \f$, i.e. the first two and the last
   * two sit on the same centre.
   *
   * The resulting rank 4 tensor is written in linearised fashion into the output
   * vector, such that later indices run faster (C-like convention).
   */
  virtual void contract_with(const iface_multivector_type& c_wb,
                             const iface_multivector_type& c_xb,
                             const iface_multivector_type& c_yb,
                             const iface_multivector_type& c_zb,
                             std::vector<Scalar>& out) const = 0;

  /** Extract a block of the untransformed tensor, i.e. in terms
   *  of the basis used by the IntegralCollection
   */
  virtual void extract_block(std::array<krims::Range<size_t>, 4> block,
                             std::vector<Scalar>& out) const = 0;

  virtual ~EriTensor_i() = default;
  EriTensor_i() = default;
  EriTensor_i(const EriTensor_i&) = default;
  EriTensor_i(EriTensor_i&&) = default;
  EriTensor_i& operator=(EriTensor_i&&) = default;
  EriTensor_i& operator=(const EriTensor_i&) = default;
};

}  // namespace gint
