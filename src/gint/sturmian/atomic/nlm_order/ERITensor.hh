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
#include "SturmintSystem.hh"
#include "gint/ERITensor_i.hh"
#include "sturmint/atomic/cs/RepulsionCalculator_i.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

template <typename RepulsionCalculator>
class ERITensor final : public ERITensor_i<scalar_type> {
 public:
  /** The precision used at the level of the sturmint library */
  typedef sturmint::real_type working_scalar_type;

  static_assert(
        std::is_base_of<sturmint::RepulsionCalculator_i, RepulsionCalculator>::value,
        "Replusion calculator needs to be derived off sturmint::RepulsionCalculator_i");

  ERITensor(const RepulsionCalculator& calculator, const SturmintSystem& system)
        : m_system_ptr("ERITensor", system), m_calculator_ptr("ERITensor", calculator) {}

  size_t n_bas() const override { return m_system_ptr->n_bas(); }

 protected:
  void compute_kernel(const std::array<krims::Range<size_t>, 4>& block,
                      kernel_type kernel) const override;

 private:
  krims::SubscriptionPointer<const SturmintSystem> m_system_ptr;
  krims::SubscriptionPointer<const RepulsionCalculator> m_calculator_ptr;
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
