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
#include <array>
#include <linalgwrap/MultiVector.hh>

namespace gint {

/** Interface to the electron repulsion integral tensor.
 *
 * Throughout this class we use chemists/Mullikan notation for the
 * integrals, i.e. element (a,b,c,d) is the integral (a b | c d)
 * with a and b being basis functions of electron 1 and c and being
 * basis functions of electron 2.
 *
 * For the tensor elements in linearised storage a C-like convention
 * is employed, i.e. later indices run faster.
 */
template <typename Scalar>
class ERITensor_i {
 public:
  typedef const linalgwrap::MutableMemoryVector_i<Scalar> iface_vector_type;
  typedef linalgwrap::MultiVector<iface_vector_type> iface_multivector_type;

  //@{
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

  virtual void contract_with(const iface_multivector_type& c_wb,
                             const iface_multivector_type& c_xb,
                             const iface_multivector_type& c_yb,
                             const iface_multivector_type& c_zb,
                             std::vector<Scalar>& out) const;
  //@}

  /** Extract a block of the untransformed tensor, i.e. in terms
   *  of the basis used by the IntegralCollection
   */
  virtual void extract_block(const std::array<krims::Range<size_t>, 4>& block,
                             std::vector<Scalar>& out) const;

  virtual ~ERITensor_i() = default;
  ERITensor_i() = default;
  ERITensor_i(const ERITensor_i&) = default;
  ERITensor_i(ERITensor_i&&) = default;
  ERITensor_i& operator=(ERITensor_i&&) = default;
  ERITensor_i& operator=(const ERITensor_i&) = default;

  /** Return the number of basis functions */
  virtual size_t n_bas() const = 0;

 protected:
  /** The type of the kernel functions, i.e a function taking the batch which is to be
   *  consumed in the current call as well as a pointer into the values.
   *
   *  The indices contained in the array of krims ranges are the absolute indices
   *  and not the indices into the data pointer.
   *
   *  The values are ordered in a C-like convention, i.e. the later indices run
   *  faster. The values pointer is managed by the kernel and should not be
   *  deallocated by the called function
   */
  typedef std::function<void(const std::array<krims::Range<size_t>, 4>&, const Scalar*)>
        kernel_type;

  // TODO Note that std::function objects are slower than lambdas or explicit
  //      templated Functors, so the batch size should not be too small.

  /** Compute a particular kernel for a block of the repulsion tensor.
   *
   * The implementing method should batch the block as it likes and for each batch
   * call the kernel function. The values should be made available by a simple
   * pointer
   */
  virtual void compute_kernel(const std::array<krims::Range<size_t>, 4>& block,
                              kernel_type kernel) const = 0;

  /** Compute a kernel for the full range of the ERI tensor */
  virtual void compute_kernel(kernel_type kernel) const {
    using krims::range;
    compute_kernel({{range(n_bas()), range(n_bas()), range(n_bas()), range(n_bas())}},
                   std::move(kernel));
  }
};

}  // namespace gint
