#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the gint authors
##
## This file is part of gint.
##
## gint is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gint is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with gint. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

. update_from_sisters.lib.sh || exit 1

update_file "krims"    "external/get.lib.sh" || exit 1
update_file "krims"    "external/get_rapidcheck.sh" || exit 1
update_file "lazyten"  "external/get_krims.sh" || exit 1
update_file "gscf"     "external/get_lazyten.sh" || exit 1

update_file "krims"    "cmake/findRapidcheck.cmake" || exit 1
update_file "krims"    "cmake/findCatch.cmake" || exit 1
update_file "krims"    "cmake/IncludeKrimsCmakeModule.cmake" || exit 1

update_file "krims"    "doc/Doxyfile.in" || exit 1

update_file "lazyten"  "templates/cc.template" "keep_header" || exit 1
update_file "lazyten"  "templates/py.template" "keep_header" || exit 1
update_file "lazyten"  "templates/cmake.template" "keep_header" || exit 1
update_file "lazyten"  "templates/hh.template" "keep_header" || exit 1
update_file "lazyten"  "templates/README.md" || exit 1

update_file "krims"    ".clang-format" || exit 1
update_file "krims"    ".clang-tidy" || exit 1
update_file "lazyten"  "update_from_sisters.lib.sh" || exit 1
