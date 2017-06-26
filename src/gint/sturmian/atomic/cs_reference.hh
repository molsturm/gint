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
#include "gint/config.hh"

#ifdef GINT_HAVE_STURMINT
#include "nlm_order.hh"
#include <sturmint/atomic/cs_reference/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_reference {
using namespace nlm_order;
using sturmint::atomic::cs_reference::Atomic;

// This integral class uses (n,l,m)-ordering: {{n,1,nmax},{l,0,n-1},{m,-l,l}}

class IntegralCollection final : public IntegralCollectionBase {
 public:
  const static std::string id;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - structure (gint::Structure): The structure of the system (i.e. an atom)
   *   - n_max (int): The maximal principle quantum number
   *   - l_max (int): Maximal azimuthal quantum number
   *   - m_max (int): Maximal magnetic quantum number
   *   - nlmbasis (NlmCollection): The precise basis triples (nlm) to use
   */
  IntegralCollection(const krims::GenMap& parameters);

  using base_type::lookup_integral;
  integral_matrix_type lookup_integral(IntegralType type) const override;

  const std::string& basis_id() const override { return id; }
  std::string basis_name() const override {
    return "Reference implementation of atomic Coulomb Sturmians";
  }

  const ERITensor_i<scalar_type>& eri_tensor() const override {
    return *m_eri_tensor_ptr;
  }

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::unique_ptr<base_type> create(const krims::GenMap& parameters) {
    return krims::make_unique<IntegralCollection>(parameters);
  }

 private:
  /** The integral calculator object */
  Atomic m_integral_calculator;

  /** Pointer to the repulsion tensor object */
  std::unique_ptr<ERITensor_i<scalar_type>> m_eri_tensor_ptr;
};

}  // namespace cs_reference
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
