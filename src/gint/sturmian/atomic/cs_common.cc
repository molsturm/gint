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

#include "cs_common.hh"

namespace gint {
namespace sturmian {
namespace atomic {

nlmCollection::nlmCollection(int nmax, int lmax, int mmax)
      : nmax(nmax), lmax(lmax), mmax(mmax) {
  assert_greater(0, nmax);
  assert_range(0, lmax, nmax);
  assert_range(0, mmax, lmax + 1);

  for (int n = 1; n <= nmax; ++n) {
    for (int l = 0; l <= min(lmax, n - 1); ++l) {
      const int ms = min(mmax, l);
      for (int m = -ms; m <= ms; ++m) {
        this->push_back(nlm_t{n, l, m});
      }
    }
  }
}

}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
