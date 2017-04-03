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
#include "IntegralCoreBase.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

class ERICoreBase : public IntegralCoreBase {
 public:
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap                                     */
  virtual void update(const krims::GenMap& map) final override;

  const coefficients_type& coeff_bo() const { return *coefficients_occupied_ptr; }

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  ERICoreBase(const SturmintSystem& system, IntegralIdentifier id)
        : IntegralCoreBase(system, id) {}
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
