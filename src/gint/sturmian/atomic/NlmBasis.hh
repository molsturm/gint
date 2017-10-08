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
#include <krims/Subscribable.hh>
#include <limits>
#include <ostream>
#include <vector>

namespace gint {
namespace sturmian {
namespace atomic {

/** Struct representing an n,l,m triple */
struct Nlm {
  int n;
  int l;
  int m;
};

//@{
/** Comparators for Nlm objects */
inline bool operator==(const Nlm& lhs, const Nlm& rhs) {
  return lhs.n == rhs.n && lhs.l == rhs.l && lhs.m == rhs.m;
}

inline bool operator!=(const Nlm& lhs, const Nlm& rhs) {
  return lhs.n != rhs.n || lhs.l != rhs.l || lhs.m != rhs.m;
}
//@}

/** Print an Nlm triple to an output stream */
std::ostream& operator<<(std::ostream& o, const Nlm& nlm);

/** Structure representing a collection of n,l,m triples */
struct NlmBasis : public std::vector<Nlm>, public krims::Subscribable {
  static constexpr int all = std::numeric_limits<int>::max();

  /** Construct an nlm basis from a maximal value for the quantum
   * numbers n (principle quantum number), l (azimuthal quantum number)
   * and m (magnetic quantum number)
   *
   * The values are automatically truncated to the allowed range.
   * I.e. if lmax == nmax, than l will only run until lmax-1.
   **/
  NlmBasis(int n_max, int l_max = all, int m_max = all) {
    for (int n = 1; n <= n_max; ++n) add_shell(n, l_max, m_max);
  }

  /** Construct an empty collection */
  NlmBasis() : m_n_max_shell(0) {}

  /** Add a shell to the basis, i.e. start the next
   *  principle quantum number and go up to a certain
   *  lmax and mmax
   */
  void add_shell(int l_max = all, int m_max = all) {
    add_shell(m_n_max_shell + 1, l_max, m_max);
    ++m_n_max_shell;
  }

  /** Add a shell to the basis and specify the principle
   *  quantum number to use */
  void add_shell(int n, int l_max, int m_max);

  /** Determine the maximal value for n in the full basis */
  int n_max() const;

  /** Determine the maximal absolute value for m in the full basis */
  int m_max() const;

  /** Determine the maximal value for l in the full basis */
  int l_max() const;

 private:
  int m_n_max_shell;
};

/** Print an NlmBasis structure to an output stream */
std::ostream& operator<<(std::ostream& o, const NlmBasis& basis);

//
// ----------------------------------------------------------------
//

inline std::ostream& operator<<(std::ostream& o, const Nlm& nlm) {
  o << '(' << nlm.n << ',' << nlm.l << ',' << nlm.m << ')';
  return o;
}

}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
