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

# Setup optional dependencies and features
# alters these things
#
#       GINT_DEPENDENCIES		everyone needs these libraries
#       GINT_DEPENDENCIES_DEBUG		debug mode needs these extras
#       GINT_DEPENDENCIES_RELEASE	release mode needs these extras
#       GINT_DEPENDENCIES_TEST		tests need these extra libraries
#

# Module to help manage optional features
include_krims_cmake_module(ProjectFeatures)

####################
#-- C++ standard --#
####################
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 14)
	message(STATUS "Detected C++14 support: Setting GINT_HAVE_CXX14")
	set(GINT_HAVE_CXX14 ON)
endif()
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 17)
	message(STATUS "Detected C++17 support: Setting GINT_HAVE_CXX17")
	set(GINT_HAVE_CXX17 ON)
endif()

################
#-- sturmint --#
################
disable_feature(sturmint)
option(GINT_ENABLE_STURMINT "Enable the sturmint library to compute Sturmian integrals." OFF)

if (GINT_ENABLE_STURMINT)
	include_krims_cmake_module(FindPackageAutocheckoutFallback)
	find_package_autocheckout_fallback(sturmint 0.0.0)

	foreach (build ${DRB_BUILD_TYPES})
		set(GINT_DEPENDENCIES_${build} ${GINT_DEPENDENCIES_${build}} ${sturmint_${build}_TARGET})
	endforeach()

	enable_feature(sturmint)
endif()

################
#--  Libint2 --#
################
disable_feature(libint)
option(GINT_ENABLE_LIBINT "Enable the libint library to compute Gaussian integrals." OFF)
option(GINT_LIBINT_USE_SYSTEM "Enable the use of a system-provided libint library" OFF)
set(GINT_LIBINT_MAX_AM 6 CACHE STRING
"Maximal angular momentum libint can perform integrals over. \
Choose a smaller value to get a faster build."
)   # Note: ORCA uses a value of 7 (up to K) in the setting above.

if (GINT_ENABLE_LIBINT)
	# Check options:
	if (GINT_LIBINT_MAX_AM LESS 2)
		message(FATAL_ERROR "GINT_LIBINT_MAX_AM needs to be at least 2")
	endif()

	# Forward parameters to included module
	set(LIBINT_VERSION 2.3.1) # We need at least this version
	set(LIBINT_SEARCH_SYSTEM ${GINT_LIBINT_USE_SYSTEM})
	set(LIBINT_MAX_AM ${GINT_LIBINT_MAX_AM})

	include(cmake/findLibint.cmake)
	unset(LIBINT_VERSION)
	unset(LIBINT_SEARCH_SYSTEM)

	set(GINT_DEPENDENCIES ${GINT_DEPENDENCIES} ${LIBINT_TARGET})
	enable_feature(libint)
endif()

################
#--  Libcint --#
################
disable_feature(libcint)
option(GINT_ENABLE_LIBCINT "Enable the libcint library to compute Gaussian integrals." OFF)

if (GINT_ENABLE_LIBCINT)
	# Forward parameters to included module
	#set(LIBCINT_VERSION 2.8.16) # We need at least this version
	message(WARNING "Due to an upstream bug, we need to track libcint master")

	include(cmake/findLibcint.cmake)
	unset(LIBINT_VERSION)

	set(GINT_DEPENDENCIES ${GINT_DEPENDENCIES} ${LIBCINT_TARGET})
	enable_feature(libcint)
endif()

##########################
#--  Static integrals  --#
##########################
disable_feature(static_integrals)
option(GINT_ENABLE_STATIC_INTEGRALS "Enable a basis types which consist entirely of pre-computed integral data. \
(Increases binary size, but useful for testing)" OFF)

if(GINT_ENABLE_STATIC_INTEGRALS)
	message(STATUS "Enabled pre-computed static integrals")
	enable_feature(static_integrals)
endif()
