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

# Finds and sets up libcint under the target name stored in
#      LIBCINT_TARGET
# such that just linking against it as a dependency does everything
# automatically.
#
# In case the libcint library is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, libcint is automatically checked out and built.
#
# Variables which change the behaviour of this module:
#
#    LIBCINT_VERSION         The libint version to be seeked
#

#
# Options and cache variables
#
option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)
message(WARNING "AUTOCHECKOUT_FORCED should be worked into this.")

#
# -------
#

# Set cmake policy for ExternalProject module and include it
if (CMAKE_VERSION VERSION_GREATER 3.1)
	cmake_policy(SET CMP0054 NEW)
endif()
include(ExternalProject)

function(SETUP_LIBCINT_FOR_EXTERNAL_BUILD TARGET)
	if ("${LIBCINT_VERSION}" STREQUAL "")
		set(LIBCINT_TAG "master")
	else()
		set(LIBCINT_TAG "v${LIBCINT_VERSION}")
	endif()

	# Special compiler flags:
	set(LIBCINT_C_FLAGS "${CMAKE_C_FLAGS}")
	enable_if_cc_compiles(LIBCINT_C_FLAGS "-Wno-unused-variable")
	enable_if_cc_compiles(LIBCINT_C_FLAGS "-Wno-unused-function")
	enable_if_cc_compiles(LIBCINT_C_FLAGS "-Wno-extra-semi")

	ExternalProject_Add(libcint
		PREFIX "${PROJECT_BINARY_DIR}/external/libcint"
		GIT_REPOSITORY "https://github.com/sunqm/libcint"
		GIT_TAG "${LIBCINT_TAG}"
		#
		CMAKE_ARGS
			# Setup compiler and compiler options
			-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
			-DCMAKE_C_FLAGS=${LIBCINT_C_FLAGS}
			-DCMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}
			-DCMAKE_BUILD_TYPE=RELEASE
			#
			# Set the place to install the library in the end
			-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/libcint
			#
			# Build but static and with PIC enabled
			# This way we can link a single gint library shared object
			# with all the code from the dependent libraries, too.
			-DENABLE_STATIC=ON -DBUILD_SHARED_LIBS=OFF
			-DCMAKE_POSITION_INDEPENDENT_CODE=ON
			-DENABLE_TESTS=ON
	)
	# TODO Right now libcint insists on building itself under the CMAKE_BUILD_TYPE
	#      RELWITHDEBINFO. Even with an explicit -DCMAKE_BUILD_TYPE=RELEASE this
	#      cannot be changed.

	ExternalProject_Get_Property(libcint install_dir)
	set(LINKFILE "${install_dir}/lib/libcint.a")

	# Setup target as an external imported library and set the include dir.
	add_library(${TARGET} STATIC IMPORTED GLOBAL)
	set_target_properties(${TARGET} PROPERTIES
		IMPORTED_LOCATION "${LINKFILE}"
	)
	include_directories(SYSTEM "${install_dir}/include")

	# Make sure that linking to ${TARGET} causes
	# a build of the libint library as well
	add_dependencies(${TARGET} libcint)
endfunction()


if (TARGET "${LIBCINT_TARGET}")
	message(STATUS "Found target ${LIBCINT_TARGET}, assume libcint already configured for build.")
	return()
endif()


if(AUTOCHECKOUT_MISSING_REPOS)
	set(LIBCINT_TARGET "External::libcint")
	SETUP_LIBCINT_FOR_EXTERNAL_BUILD(${LIBCINT_TARGET})

	message(STATUS "Setting up libcint as an external project")
	return()
else()
	message(FATAL_ERROR "For libcint currently AUTOCHECKOUT_MISSING_REPOS is required.")
	return()
endif()

