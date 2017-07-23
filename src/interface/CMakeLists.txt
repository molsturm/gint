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

#
# Parse build type
#
set(GINT_INTERFACE_BUILD_TYPE "Release" CACHE STRING
"The build type to use for the gint bindings to python and other languages. \
Valid are Debug and Release.")

set(module_BUILD_TYPE "DEBUG")
if (GINT_INTERFACE_BUILD_TYPE STREQUAL "Release")
	set(interfaces_BUILD_TYPE "RELEASE")
elseif(GINT_INTERFACE_BUILD_TYPE  STREQUAL "Debug")
	set(interfaces_BUILD_TYPE "DEBUG")
else()
	message(FATAL_ERROR
		"The value of GINT_INTERFACE_BUILD_TYPE needs to be either \"Debug\" or \"Release\".")
endif()

#
# Pass onto subdirectories
#
if(GINT_INTERFACE_PYTHON)
	add_subdirectory(python)
endif()
