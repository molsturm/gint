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

## DO NOT EDIT
## This file is automatically generated from a file in the repository "lazyten".
## Edit the original and call the script "update_from_sister_repos.sh" instead.

# Finds and sets up krims under the target names stored in
#      krims_DEBUG_TARGET     (Debug version)
#      krims_RELEASE_TARGET   (Release version)
# such that just linking against it as a dependency does everything
# automatically.
#
# In case the krims library is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, krims is automatically checked out and built.
# Otherwise a fatal error is produced.
#

#
# Options and properties required
#
option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)

#
# -------
#

if (TARGET "${krims_DEBUG_TARGET}"  OR TARGET "${krims_RELEASE_TARGET}")
	message(STATUS "Found krims targets, assume krims already configured for build.")
	return()
endif()

# Try to find krims somewhere
find_package(krims ${KRIMS_VERSION} QUIET CONFIG)
mark_as_advanced(krims_DIR)

if ("${krims_DIR}" STREQUAL "krims_DIR-NOTFOUND")
	if (AUTOCHECKOUT_MISSING_REPOS)
		execute_process(
			COMMAND "sh" "get_krims.sh"
			WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/external"
			RESULT_VARIABLE RES
		)
		if (NOT RES EQUAL 0)
			message(FATAL_ERROR "Getting krims from git failed with error: ${RES}")
		endif()

		#
		# Proceed to configure krims
		#
		add_subdirectory(${PROJECT_SOURCE_DIR}/external/krims)
		include_directories(${PROJECT_SOURCE_DIR}/external/krims/src)
		include_directories(${PROJECT_BINARY_DIR}/external/krims/src)

		# Extract version from CMakeLists.txt:
		file(STRINGS "${PROJECT_SOURCE_DIR}/external/krims/CMakeLists.txt"
			VERSION_RAW
			REGEX "krims VERSION [0-9.]+"
			LIMIT_COUNT 1)
		string(REGEX MATCH "[0-9.]+" GOT_VERSION "${VERSION_RAW}")

		# Compare against what is needed
		if("${GOT_VERSION}" VERSION_LESS "${KRIMS_VERSION}")
			message(FATAL_ERROR "Inconsistency in the repo: \
Version ${KRIMS_VERSION} of krims was requested, but only version ${GOT_VERSION} \
was found.")
		endif()

		return()
	endif()

	message(FATAL_ERROR "Could not find krims library.
Either provide the installation prefix of the krims library in the environment \
variable krims_DIR or enable autocheckout via '-DAUTOCHECKOUT_MISSING_REPOS=ON'.")
endif()

message(WARNING "This part of findKrims has never been tested.")

# Setup library targets
set(krims_DEBUG_TARGET   "Upstream::krims.g"
	CACHE INTERNAL "Target name of debug version of krims")
set(krims_RELEASE_TARGET "Upstream::krims"
	CACHE INTERNAL "Target name of release version of krims")

# Check that all required targets are available.
foreach(build ${DRB_BUILD_TYPES})
	if(NOT TARGET "${krims_${build}_TARGET}")
		message(FATAL_ERROR "We could not find a ${build} version of krims at this location. \
		Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
		rebuild krims with a ${build} version as well.")
	endif()
endforeach()

#TODO check that we don't need extra stuff like in findLinalgwrap in gscf

message(STATUS "Found krims config at ${krims_CONFIG}")
