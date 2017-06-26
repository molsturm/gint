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

#include "SturmintSystem.hh"
#include <complex>
#ifdef GINT_HAVE_STURMINT

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

SturmintSystem::SturmintSystem(scalar_type Z, scalar_type k, const NlmBasis& basis_)
      : k(std::move(k)), basis(), Z(std::move(Z)) {
  basis.reserve(basis_.size());
  for (const Nlm& nlm : basis_) {
    assert_greater(0, nlm.n);
    assert_greater_equal(0, nlm.l);
    assert_greater_equal(std::abs(nlm.m), nlm.l);
    basis.push_back(nlm_type{nlm.n, nlm.l, nlm.m});
  }
}

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
