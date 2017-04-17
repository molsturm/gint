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

#include "qnum.hh"
#include <algorithm>
#include <locale>

namespace gint {
namespace qnum {

namespace {
const std::array<char, 13> map = {
      {'S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'}};
}

int letter_to_am(char l) {
  auto it = std::find(std::begin(map), std::end(map), std::toupper(l));
  assert_throw(it != std::end(map), ExcInvalidAmLetter(l));
  return static_cast<int>(it - std::begin(map));
}

char am_to_letter(int l) {
  assert_range(static_cast<long>(0), static_cast<long>(l), static_cast<long>(map.size()));
  return map.at(static_cast<size_t>(l));
}

}  // namespace qnum
}  // namespace gint
