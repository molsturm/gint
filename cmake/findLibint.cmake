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
# Options and properties required
#
option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)

#
# -------
#

function(SETUP_LIBINT2_FOR_EXTERNAL_BUILD TARGET)
	# Determine compiler flags which are in use in outer project
	# and remove all -W and -f flags
	string(REGEX REPLACE "(-(W|f)[^ ]+|-pedantic)" "" TMP
		${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS})
	string(REGEX REPLACE "  *"  " "  INNER_COMPILE_OPTS ${TMP})

	set(CONFIGURE_OPTS
		# Optimisation flags for the inner compiler
		#      => Forward what we have in the global project.
		"--with-cxxgen-optflags=${INNER_COMPILE_OPTS} -std=c++${CMAKE_CXX_STANDARD}"
	)

	if (BUILD_SHARED_LIBS)
		# Enable position-independent code such that we can link libint
		# into a shared libgint.so object.
		set (CONFIGURE_OPTS "--with-pic" ${CONFIGURE_OPTS})
	endif()

	include(DefaultExternalProjects)
	setup_autotools_project(libint
		CONFIGURE_OPTS ${CONFIGURE_OPTS}
		GIT_REPOSITORY "https://github.com/evaleev/libint.git"
		GIT_TAG "v2.3.0-beta.3"
		#
		# TODO Find out how to test libint ... e.g. TEST_COMMAND make check

		CHECKPOINT_TO "<INSTALL_DIR>/lib/libint2.a"
		CHECKPOINT_FROM "<BINARY_DIR>/lib/.libs/libint2.a"
		# This has to be the last option!
		CONFIGURE_OPTS ${CONFIGURE_OPTS}
	)
	ExternalProject_Get_Property(libint install_dir)
	set(LINKFILE "${install_dir}/lib/libint2.a")

	# Setup target as an external imported library and set the include dir.
	add_library(${TARGET} INTERFACE IMPORTED GLOBAL)
	set_target_properties(${TARGET} PROPERTIES
		INTERFACE_LINK_LIBRARIES "${LINKFILE}"
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

	# Setup link target
	add_library(${TARGET} INTERFACE IMPORTED GLOBAL)
	set_target_properties(${TARGET} PROPERTIES
		INTERFACE_LINK_LIBRARIES ${LIBINT_LIBRARIES}
		INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${LIBINT_INCLUDE_DIRS}
	# TODO properly adapt this!
	#	INTERFACE_COMPILE_DEFINITIONS "-DDATADIR=../../../../../libint_install/share/libint/2.2.0-beta1/basis/"
	)

	set(LIBINT_VERSION ${LIBINT_VERSION} PARENT_SCOPE)
endfunction()

#
# -------
#

if (TARGET "${LIBINT_TARGET}")
	message(STATUS "Found libint target, assume libint already confgured for build.")
	return()
endif()

# Since libint2 needs the Eigen matrix library, we try to find that first:
# Note: Eigen is available under Eigen3::Eigen
find_package(Eigen3 REQUIRED)

if (LIBINT_SEARCH_SYSTEM STREQUAL "" OR LIBINT_SEARCH_SYSTEM)
	set(LIBINT_TARGET "System::libint")
	SETUP_SYSTEM_LIBINT(${LIBINT_TARGET} ${LIBINT_VERSION})
	if (LIBINT_FOUND)
		target_link_libraries(${LIBINT_TARGET} INTERFACE Eigen3::Eigen)

		message(STATUS "Found libint version ${LIBINT_VERSION} at ${LIBINT_LIBRARY}")
		return()
	endif()
endif()
if(AUTOCHECKOUT_MISSING_REPOS)
	set(LIBINT_TARGET "External::libint")
	SETUP_LIBINT2_FOR_EXTERNAL_BUILD(${LIBINT_TARGET})
	target_link_libraries(${LIBINT_TARGET} INTERFACE Eigen3::Eigen)

	message(STATUS "Setting up libint2 as an external project")
	return()
endif()

message(FATAL_ERROR "Could not find libint library.
Either provide hints to the installation using the cmake variables LIBINT_INCLUDE_DIR
and LIBINT_LIBRARY or enable autocheckout via '-DAUTOCHECKOUT_MISSING_REPOS=ON'.")
