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
#include <gint/IntegralCoreBase.hh>
#include <gint/config.hh>
#include <krims/GenMap.hh>
#include <linalgwrap/MultiVector.hh>

namespace gint {

template <typename StoredMatrix>
class CoefficientContainer {
 public:
  typedef typename StoredMatrix::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  /** \brief Update the internal data
   *
   * Listens to the keys provided in IntegralUpdeteKeys.hh
   */
  void update(const krims::GenMap& map);

  //! Current occupied coefficients as a reference
  const coefficients_type& coeff_bo() const {
    assert_dbg(coeff_bo_ptr != nullptr, krims::ExcInvalidPointer());
    return *coeff_bo_ptr;
  }

  //! Current occupied coefficients as a pointer
  coefficients_ptr_type coeff_bo_ptr;

 protected:
  // TODO The better solution for this would be a lazy matrix for the density object.
  real_valued::stored_matrix_type compute_density_matrix() const {
    auto dens = linalgwrap::outer_prod_sum(coeff_bo(), coeff_bo());
    assert_dbg(dens.is_symmetric(), krims::ExcInternalError());
    return dens;
  }
};

}  // namespace gint
