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

#include "ERICore.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

/** A less-optimised version of ERICore, which uses the working_scalar_type
 *  as defined by the sturmint real_type in order to perform apply
 *  and other operations as well.
 *
 *  The class is only instantiated for some RepulsionCalculator classes.
 *  See the cc file for details
 */
template <typename RepulsionCalculator>
class ERICoreHighPrecision : public ERICore<RepulsionCalculator> {
 public:
  using base_type = ERICore<RepulsionCalculator>;
  using base_type::ERICore;
  using working_scalar_type = typename base_type::working_scalar_type;
  using base_core_type      = typename base_type::base_core_type;

  /** Apply to a multivector, internally use working_scalar_type to do computation */
  void apply(const const_multivector_type& x, multivector_type& y,
             const lazyten::Transposed mode = lazyten::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** Clone the matrix expression */
  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new ERICoreHighPrecision(*this));
  }
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
