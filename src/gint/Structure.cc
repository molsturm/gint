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

#include "Structure.hh"

namespace gint {

std::ostream& operator<<(std::ostream& o, const Structure& st) {
  auto it = std::begin(st);
  if (it != std::end(st)) o << *(it++);
  for (; it != std::end(st); ++it) o << '\n' << *it;
  return o;
}

}  // namespace gint
