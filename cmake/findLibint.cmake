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

# Finds and sets up libint under the target name stored in
#      LIBINT_TARGET
# such that just linking against it as a dependency does everything
# automatically.
#
# In case the libint library is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, libint is automatically checked out and built.
#
# Variables which change the behaviour of this module:
#
#    LIBINT_VERSION         The libint version to be seeked
#    LIBINT_SEARCH_SYSTEM   Allow to search the system for an installed libint.
#    LIBINT_MAX_AM          Maximal angular momentum to support inside libint.

#
# ----------------------------------------------------------------------
#

#
# Options and cache variables
#
option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)
message(WARNING "AUTOCHECKOUT_FORCED should be worked into this.")

#
# -------
#

function(SETUP_LIBINT2_FOR_EXTERNAL_BUILD TARGET LIBINT_MAX_AM)
	# Determine compiler flags which are in use in outer project
	# and remove all -W and -f flags
	string(REGEX REPLACE "(-(W|f)[^ ]+|-pedantic)" "" TMP
		${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS})
	string(REGEX REPLACE "  *"  " "  INNER_COMPILE_OPTS ${TMP})

	set(LIBINT_OPT_AM 4)
	if (LIBINT_MAX_AM LESS 4)
		set(LIBINT_OPT_AM ${LIBINT_MAX_AM})
	endif()

	if ("${LIBINT_VERSION}" STREQUAL "")
		set(LIBINT_TAG "master")
	else()
		set(LIBINT_TAG "v${LIBINT_VERSION}")
	endif()

	# Get the flag for the C++ standard
	cxx_standard_flag(${CMAKE_CXX_STANDARD} CXX_STANDARD_FLAG)

	# Get flags for the standard library to use
	stdlib_cxx_flags(${DRB_CXX_STANDARD_LIBRARY} CLANG_CXX_FLAGS CLANG_LD_FLAGS)

	set(CONFIGURE_OPTS
		# Export compiler selection
		#
		"--with-ranlib=${CMAKE_RANLIB}"
		"--with-ar=${CMAKE_AR}"
		"--with-ld=${CMAKE_LD}"
		"CC=${CMAKE_C_COMPILER}"
		"CXX=${CMAKE_CXX_COMPILER}"
		#
		# Optimisation flags for all compilation processes inside libint
		"CXXFLAGS=${CXX_STANDARD_FLAG} ${CLANG_CXX_FLAGS}"
		"LDFLAGS=${CLANG_LD_FLAGS}"
		#
		# Optimisation flags for the inner compiler
		#      => Forward what we have in the global project.
		"--with-cxxgen-optflags=${INNER_COMPILE_OPTS}"
		#
		# Enable integrals up to the maximal angular momentum requested
		"--with-max-am=${LIBINT_MAX_AM}"
		#
		# Optimise maximally up to angular momentum 4
		"--with-opt-am=${LIBINT_OPT_AM}"
		#
		# Disable unrolling shell sets into integrals
		"--disable-unrolling"
		#
		# Enable hand-written generic code for other integrals
		"--enable-generic-code"
	)

	# TODO I do not quite understand why, but we seem to need this regardles
	#      of whether we end up linking statically or dynamically for the
	#      test executables to link without a relocation error
	#if (BUILD_SHARED_LIBS)
		# Enable position-independent code such that we can link libint
		# into a shared libgint.so object.
		set (CONFIGURE_OPTS "--with-pic" ${CONFIGURE_OPTS})
	#endif()

	include(DefaultExternalProjects)
	setup_autotools_project(libint
		GIT_REPOSITORY "https://github.com/evaleev/libint.git"
		GIT_TAG "${LIBINT_TAG}"
		#
		# TODO Find out how to test libint ... e.g. TEST_COMMAND make check

		# This only works once cmake >= 3.3:
		#CHECKPOINT_TO "<INSTALL_DIR>/lib/libint2.a"
		CHECKPOINT_TO "${PROJECT_BINARY_DIR}/external/libint/lib/libint2.a"
		CHECKPOINT_FROM "<BINARY_DIR>/lib/.libs/libint2.a"
		# This has to be the last option!
		CONFIGURE_OPTS ${CONFIGURE_OPTS}
	)
	ExternalProject_Get_Property(libint install_dir)
	set(LINKFILE "${install_dir}/lib/libint2.a")

	# Setup target as an external imported library and set the include dir.
	add_library(${TARGET} STATIC IMPORTED GLOBAL)
	set_target_properties(${TARGET} PROPERTIES
		IMPORTED_LOCATION "${LINKFILE}"
	)
	include_directories(SYSTEM "${install_dir}/include")

	# Make sure that linking to ${TARGET} causes
	# a build of the libint library as well
	add_dependencies(${TARGET} libint)
endfunction()

function(SETUP_SYSTEM_LIBINT TARGET VERSION)
	find_package(Libint2 ${VERSION} QUIET MODULE)
	if (NOT LIBINT_FOUND)
		return()
	endif()

	# Run find_package again, but this time do not be quiet:
	find_package(Libint2 ${VERSION} REQUIRED MODULE)

	# TODO Determine max angular momentum and check

	# Setup link target
	add_library(${TARGET} INTERFACE IMPORTED GLOBAL)
	set_target_properties(${TARGET} PROPERTIES
		INTERFACE_LINK_LIBRARIES "${LIBINT_LIBRARIES}"
		INTERFACE_INCLUDE_DIRECTORIES "${LIBINT_INCLUDE_DIRS}"
	)

	set(LIBINT_VERSION ${LIBINT_VERSION} PARENT_SCOPE)
	set(LIBINT_FOUND ${LIBINT_FOUND} PARENT_SCOPE)
endfunction()

#
# -------
#

if (TARGET "${LIBINT_TARGET}")
	message(STATUS "Found libint target, assume libint already confgured for build.")
	return()
endif()

# Since libint2 needs the Eigen matrix library, we try to find that first:
# Note: Eigen is exposed as Eigen3::Eigen by find_package
set(EIGEN_TARGET "Eigen3::Eigen")
find_package(Eigen3 QUIET)
if (${Eigen3_DIR} MATCHES "-NOTFOUND" OR NOT TARGET ${EIGEN_TARGET})
	set(EIGEN_TARGET "Eigen")
	add_library(${EIGEN_TARGET} INTERFACE IMPORTED)
	find_path(Eigen3_DIR signature_of_eigen3_matrix_library
		PATH_SUFFIXES eigen3
		DOC "The eigen3 include directory."
	)
	if (${Eigen3_DIR} MATCHES "-NOTFOUND")
		message(FATAL_ERROR "Could not find the eigen3 library anywhere.
Try setting Eigen3_DIR to the location.")
	endif()
	set_target_properties(${EIGEN_TARGET} PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES ${Eigen3_DIR}
	)
	message(STATUS "Found eigen3 include directory at ${Eigen3_DIR}")
endif()

if (LIBINT_SEARCH_SYSTEM STREQUAL "" OR LIBINT_SEARCH_SYSTEM)
	set(LIBINT_TARGET "System::libint")
	SETUP_SYSTEM_LIBINT(${LIBINT_TARGET} ${LIBINT_VERSION})
	if (LIBINT_FOUND)
		target_link_libraries(${LIBINT_TARGET} INTERFACE ${EIGEN_TARGET})

		message(STATUS "Found libint version ${LIBINT_VERSION} at ${LIBINT_LIBRARY}")
		return()
	endif()
endif()
if(AUTOCHECKOUT_MISSING_REPOS)
	set(LIBINT_TARGET "External::libint")
	SETUP_LIBINT2_FOR_EXTERNAL_BUILD(${LIBINT_TARGET} ${LIBINT_MAX_AM})
	target_link_libraries(${LIBINT_TARGET} INTERFACE ${EIGEN_TARGET})

	message(STATUS "Setting up libint2 as an external project")
	return()
endif()

message(FATAL_ERROR "Could not find libint library. \
Either enable the use of a system-provided libint (using LIBINT_SEARCH_SYSTEM, \
provide hints to the installation using the cmake variables LIBINT_INCLUDE_DIR \
and LIBINT_LIBRARY or enable autocheckout via '-DAUTOCHECKOUT_MISSING_REPOS=ON'.")
