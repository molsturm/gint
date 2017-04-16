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
#include "gint/config.hh"
#include <krims/FileUtils/FindDataFile.hh>

namespace gint {

/** Specialisation of the krims::FindDataFile functor for the gint library */
struct FindDataFile : public krims::FindDataFile {
  FindDataFile() : krims::FindDataFile("gint") {
    extra_directories.push_back(std::string(detail::data_download_dir));
    extra_directories.push_back(std::string(detail::data_install_dir));
  }
};

/** Search a datafile in gint's static data directories
 *  and return the full path to it.
 *
 *  In case the file could not be found a krims::ExcDatafileNotFound
 *  exception is raised.
 */
inline std::string find_data_file(const std::string& file) {
  return FindDataFile{}(file);
}

}  // namespace gint
