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

#include "Basis.hh"
#include "Shell.hh"
#include <iterator>

namespace gint {
namespace gaussian {

Basis::Basis(const Structure& structure, const BasisSet& set) {
  for (const Atom& a : structure) {
    // Get shells for atom and copy them inside:
    const auto& shells = set.atomic_number_to_shells.at(a.atomic_number);

    auto set_origin = [&a](Shell s) {
      s.origin = a.coords;
      return s;
    };

    std::transform(std::begin(shells), std::end(shells),
                   std::back_insert_iterator<Basis>(*this), set_origin);
  }
}

}  // namespace gaussian
}  // namespace gint
