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
#include "NlmBasis.hh"
#include "gint/config.hh"
#include <krims/Subscribable.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

// Everything in this namespace is real-valued:
using namespace gint::real_valued;

/** Structure, which collects the system information relevant to
 *  compute sturmint integrals */
struct SturmintSystem : public krims::Subscribable {
  scalar_type k;
  scalar_type Z;
  NlmBasis basis;

  /** Return the number of basis functions */
  size_t n_bas() const { return basis.size(); }
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
