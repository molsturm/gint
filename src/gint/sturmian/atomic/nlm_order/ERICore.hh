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

#include "IntegralCoreBase.hh"
#include "gint/CoefficientContainer.hh"
#include "sturmint/atomic/cs/RepulsionCalculator_i.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

/** Integral core class for computing sturmian repulsion integrals in nlm order.
 *
 *  The class is only instantiated for some RepulsionCalculator classes.
 *  See the cc file for details
 */
template <typename RepulsionCalculator>
class ERICore : public IntegralCoreBase, public CoefficientContainer<stored_matrix_type> {
  // Note: The CoefficientContainer contains the current coefficients and allows access
  // to the density derived from these.
 public:
  /** The precision used at the level of the sturmint library */
  typedef sturmint::real_type working_scalar_type;

  static_assert(
        std::is_base_of<sturmint::RepulsionCalculator_i, RepulsionCalculator>::value,
        "Replusion calculator needs to be derived off sturmint::RepulsionCalculator_i");

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap                                     */
  virtual void update(const krims::GenMap& map) final override {
    CoefficientContainer<stored_matrix_type>::update(map);
    assert_dbg(coeff_bo_ptr == nullptr || coeff_bo().n_vectors() == 0 ||
                     coeff_bo_ptr->n_elem() == n_bas(),
               krims::ExcSizeMismatch(coeff_bo_ptr->n_elem(), n_bas()));
  }

  ERICore(const RepulsionCalculator& calculator, const SturmintSystem& system,
          IntegralIdentifier id)
        : IntegralCoreBase(system, id), m_calculator_ptr("ERICore", calculator) {
    assert_internal(id.integral_type() == IntegralType::exchange ||
                    id.integral_type() == IntegralType::coulomb);
  }

  /** \brief Compute and return an element of the matrix
     J_{ab} = J_{abcd} Cocc_{cp} Cocc_{dp} = J_{abcd} P_{cd}
     K_{ab} = J_{cbad} Cocc_{cp} Cocc_{dp} = J_{acbd} P_{cd}
  */
  scalar_type operator()(size_t a, size_t b) const final override {
    return static_cast<scalar_type>(compute_jk_element(a, b));
  }

  /** Clone the matrix expression */
  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new ERICore(*this));
  }

 protected:
  /** Compute a single element at the precision offered by the sturmian library */
  working_scalar_type compute_jk_element(size_t a, size_t b) const;

 private:
  krims::SubscriptionPointer<const RepulsionCalculator> m_calculator_ptr;
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
