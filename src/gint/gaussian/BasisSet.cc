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

DefaultAngularFunctions lookup_default_angular_functions(const std::string& name) {
  const std::string nrm = detail::normalise_basisset_name(name);
  // We need cartesian ds for:
  //     - 3-21g
  //     - 4-31g
  //     - 6-21g
  //     - 6-31??g??
  // BUT NOT
  //     - 6-311??g??

  const size_t gpos = nrm.find('g');
  if (gpos == std::string::npos) {
    // No g found, so not one of the Cartesian guys
    return DefaultAngularFunctions::Pure;
  }

  if (nrm.find("3-21") == 0 || nrm.find("4-31") == 0 || nrm.find("6-21") == 0 ||
      nrm.find("6-31") == 0) {
    if (nrm.size() == 4 || (nrm[4] != '+' && nrm[4] != 'g')) {
      // To make sure 6-311??g?? get pure ds
      return DefaultAngularFunctions::Pure;
    }
    if (gpos + 1 == nrm.size() || nrm[gpos + 1] == '*') {
      return DefaultAngularFunctions::CartesianForD;
    }
  }
  return DefaultAngularFunctions::Pure;
}

BasisSet lookup_basisset(const std::string& name) {
  const std::string datafile = "basis/" + detail::normalise_basisset_name(name) + ".g94";

  // Find files, but also search in current working directory (i.e. basis files in
  // the basis subfolder will be recognised)
  gint::FindDataFile find;
  find.cwd_suffixes.push_back(".");
  find.env_vars.push_back("BASIS_DIR");

  BasisSet b = read_basisset(find(datafile), BasisSetFileFormat::Gaussian94);
  b.name = detail::normalise_basisset_name(name, false);

  const DefaultAngularFunctions func = lookup_default_angular_functions(name);
  if (func == DefaultAngularFunctions::Pure) {
    // This is done by default, so
    return b;
  }

  for (auto& kv : b.atomic_number_to_shells) {
    for (Shell& s : kv.second) {
      if (func == DefaultAngularFunctions::Cartesian) {
        s.pure = false;
      } else if (func == DefaultAngularFunctions::CartesianForD) {
        s.pure = (s.l > 2);
      }
    }
  }
  return b;
}

}  // namespace gaussian
}  // namespace gint
