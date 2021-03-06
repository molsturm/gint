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

#include "NlmBasis.hh"
#include <krims/Algorithm.hh>
#include <krims/ExceptionSystem.hh>

namespace gint {
namespace sturmian {
namespace atomic {

void NlmBasis::add_shell(int n, int l_max, int m_max) {
  assert_greater(0, n);
  assert_greater_equal(0, l_max);
  assert_greater_equal(0, m_max);

  for (int l = 0; l <= std::min(l_max, n - 1); ++l) {
    const int ms = std::min(m_max, l);
    for (int m = -ms; m <= ms; ++m) {
      this->push_back({n, l, m});
    }
  }
}

int NlmBasis::n_max() const {
  auto compare_n = [](const Nlm& lhs, const Nlm& rhs) { return lhs.n < rhs.n; };
  return std::max_element(begin(), end(), std::move(compare_n))->n;
}

/** Determine the maximal absolute value for m in the full basis */
int NlmBasis::m_max() const {
  auto compare_absm = [](const Nlm& lhs, const Nlm& rhs) {
    return std::abs(lhs.m) < std::abs(rhs.m);
  };
  return std::abs(std::max_element(begin(), end(), std::move(compare_absm))->m);
}

/** Determine the maximal value for l in the full basis */
int NlmBasis::l_max() const {
  auto compare_l = [](const Nlm& lhs, const Nlm& rhs) { return lhs.l < rhs.l; };
  return std::max_element(begin(), end(), std::move(compare_l))->l;
}

std::ostream& operator<<(std::ostream& o, const NlmBasis& basis) {
  o << krims::join(std::begin(basis), std::end(basis), ", ");
  return o;
}

}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
