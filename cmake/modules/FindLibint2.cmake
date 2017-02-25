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

# Find module for libint2
#
# Populates the variables:
#     LIBINT_FOUND - set to true if the library is found
#     LIBINT_INCLUDE_DIRS - list of required include directories
#     LIBINT_LIBRARIES - list of libraries to be linked
#     LIBINT_MAJOR_VERSION - major version number
#     LIBINT_MINOR_VERSION - minor version number
#     LIBINT_MICRO_VERSION - micro version number
#     LIBINT_VERSION_STRING - version number as a string (ex: "1.0.4")
#

# Empty variables
set(LIBINT_TARGET "")
set(LIBINT_MAJOR_VERSION 0)
set(LIBINT_MINOR_VERSION 0)
set(LIBINT_MICRO_VERSION 0)
set(LIBINT_VERSION_STRING "")

# Try to find lib and include
find_library(LIBINT_LIBRARY NAMES int2)
find_path(LIBINT_INCLUDE_DIR NAMES libint2.hpp)

if(NOT LIBINT_INCLUDE_DIR MATCHES "NOTFOUND"
	AND EXISTS "${LIBINT_INCLUDE_DIR}/libint2/config.h")
	# Extract version information from config header

	file(STRINGS "${LIBINT_INCLUDE_DIR}/libint2/config.h" _header_contents
		REGEX "#define LIBINT_[A-Z]+_VERSION ")
	string(REGEX REPLACE ".*#define LIBINT_MAJOR_VERSION ([0-9]+).*" "\\1"
		LIBINT_MAJOR_VERSION "${_header_contents}")
	string(REGEX REPLACE ".*#define LIBINT_MINOR_VERSION ([0-9]+).*" "\\1"
		LIBINT_MINOR_VERSION "${_header_contents}")
	string(REGEX REPLACE ".*#define LIBINT_MICRO_VERSION ([0-9]+).*" "\\1"
		LIBINT_MICRO_VERSION "${_header_contents}")
	unset(_header_contents)

	set(LIBINT_VERSION_STRING "${LIBINT_MAJOR_VERSION}.${LIBINT_MINOR_VERSION}.${LIBINT_MICRO_VERSION}")
endif()

# Checks 'REQUIRED', 'QUIET' and versions.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libint
  REQUIRED_VARS LIBINT_LIBRARY LIBINT_INCLUDE_DIR
  VERSION_VAR LIBINT_VERSION_STRING
  FAIL_MESSAGE "Could NOT find libint2 on the system"
)

if (LIBINT_FOUND)
	set(LIBINT_INCLUDE_DIRS ${LIBINT_INCLUDE_DIR})
	set(LIBINT_LIBRARIES ${LIBINT_LIBRARY})
endif ()

# Hide internal variables
mark_as_advanced(LIBINT_INCLUDE_DIR LIBINT_LIBRARY)
