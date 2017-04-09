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
#include "gint/sturmian/atomic/NlmBasis.hh"
#include <krims/Subscribable.hh>
#include <sturmint/atomic/cs/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

// Import real-valued structures and nlm type:
using namespace gint::real_valued;
using nlm_type = sturmint::atomic::cs::nlm_t;

/** Structure, which collects the system information relevant to
 *  compute sturmint integrals */
struct SturmintSystem : public krims::Subscribable {
  scalar_type k;
  std::vector<nlm_type> basis;
  scalar_type Z;

  /** Return the number of basis functions */
  size_t n_bas() const { return basis.size(); }

  SturmintSystem() = default;
  SturmintSystem(scalar_type Z, scalar_type k, std::vector<nlm_type> basis)
        : k(std::move(k)), basis(std::move(basis)), Z(std::move(Z)) {}

  SturmintSystem(scalar_type Z, scalar_type k, const NlmBasis& basis_);
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
