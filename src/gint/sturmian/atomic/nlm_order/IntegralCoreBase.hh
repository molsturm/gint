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
#include "gint/IntegralCoreBase.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

/** IntegralCoreBase for sturmint integrals */
class IntegralCoreBase : public gint::IntegralCoreBase<stored_matrix_type> {
 public:
  using base_core_type = gint::IntegralCoreBase<stored_matrix_type>;

  bool has_transpose_operation_mode() const final override { return true; }
  size_t n_rows() const final override { return m_system_ptr->n_bas(); }
  size_t n_cols() const final override { return m_system_ptr->n_bas(); }
  size_t n_bas() const { return m_system_ptr->n_bas(); }
  IntegralIdentifier id() const final override { return m_id; }

  IntegralCoreBase(const SturmintSystem& system, IntegralIdentifier id)
        : m_system_ptr("sturmian::IntegralCoreBase", system), m_id(id) {}

 private:
  krims::SubscriptionPointer<const SturmintSystem> m_system_ptr;

 protected:
  IntegralIdentifier m_id;
  const SturmintSystem& system() const { return *m_system_ptr; }
  const IntegralType& type() const { return m_id.integral_type(); }
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
