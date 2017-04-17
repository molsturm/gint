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

#include "BasisSet.hh"
#include "Shell.hh"
#include "gint/find_data_file.hh"
#include "read_basisset.hh"
#include <iterator>

namespace gint {
namespace gaussian {

namespace detail {
std::string normalise_basisset_name(const std::string& name, bool full = true) {
  auto normalise_char = [&full](char c) {
    switch (c) {
      case '/':
        return full ? 'I' : '/';
      default:
        return static_cast<char>(::tolower(c));
    }
  };

  std::string normalised;
  normalised.reserve(name.length());
  std::transform(std::begin(name), std::end(name),
                 std::back_insert_iterator<std::string>(normalised), normalise_char);

  return normalised;
}
}  // namespace detail

BasisSet lookup_basisset(const std::string& name) {
  // TODO Augmentation and other stuff ??
  //      Shall we use an extra file like libint?

  const std::string datafile = "basis/" + detail::normalise_basisset_name(name) + ".g94";

  // Find files, but also search in current working directory (i.e. basis files in
  // the basis subfolder will be recognised)
  gint::FindDataFile find;
  find.cwd_suffixes.push_back(".");

  const std::string fullpath = find(datafile);
  std::ifstream f(fullpath);
  BasisSet b = read_basisset(f, BasisSetFileFormat::Gaussian94);
  b.name = detail::normalise_basisset_name(name, false);
  b.filename = fullpath;
  return b;
}

}  // namespace gaussian
}  // namespace gint
