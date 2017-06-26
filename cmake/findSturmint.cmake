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

# Finds and sets up sturmint under the target names stored in
#      sturmint_DEBUG_TARGET     (Debug version)
#      sturmint_RELEASE_TARGET   (Release version)
# such that just linking against it as a dependency does everything
# automatically.
#
# In case the sturmint library is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, sturmint is automatically checked out and built.
# Otherwise a fatal error is produced.
#

#
# Options and properties required
#
option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)

#
# -------
#

if (TARGET "${sturmint_DEBUG_TARGET}"  OR TARGET "${sturmint_RELEASE_TARGET}")
	message(STATUS "Found sturmint targets, assume sturmint already configured for build.")
	return()
endif()

# Try to find sturmint somewhere
find_package(sturmint ${STURMINT_VERSION} QUIET CONFIG)
mark_as_advanced(sturmint_DIR)

if ("${sturmint_DIR}" STREQUAL "sturmint_DIR-NOTFOUND")
	if (AUTOCHECKOUT_MISSING_REPOS)
		execute_process(
			COMMAND "sh" "get_sturmint.sh"
			WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/external"
			RESULT_VARIABLE RES
		)
		if (NOT RES EQUAL 0)
			message(FATAL_ERROR "Getting sturmint from git failed with error: ${RES}")
		endif()

		#
		# Proceed to configure sturmint
		#
		add_subdirectory(${PROJECT_SOURCE_DIR}/external/sturmint)
		include_directories(${PROJECT_SOURCE_DIR}/external/sturmint/src)
		include_directories(${PROJECT_BINARY_DIR}/external/sturmint/src)

		# Extract version from CMakeLists.txt:
		file(STRINGS "${PROJECT_SOURCE_DIR}/external/sturmint/CMakeLists.txt"
			VERSION_RAW
			REGEX "sturmint VERSION [0-9.]+"
			LIMIT_COUNT 1)
		string(REGEX MATCH "[0-9.]+" GOT_VERSION "${VERSION_RAW}")

		# Compare against what is needed
		if("${GOT_VERSION}" VERSION_LESS "${STURMINT_VERSION}")
			message(FATAL_ERROR "Inconsistency in the repo: \
Version ${STURMINT_VERSION} of sturmint was requested, but only version ${GOT_VERSION} \
was found.")
		endif()

		return()
	endif()

	message(FATAL_ERROR "Could not find sturmint library.
Either provide the installation prefix of the sturmint library in the environment \
variable sturmint_DIR or enable autocheckout via -DAUTOCHECKOUT_MISSING_REPOS=ON.")
endif()

message(WARNING "This part of findSturmint has never been tested.")

# Setup library targets
set(sturmint_DEBUG_TARGET   "Upstream::sturmint.g"
	CACHE INTERNAL "Target name of debug version of sturmint")
set(sturmint_RELEASE_TARGET "Upstream::sturmint"
	CACHE INTERNAL "Target name of release version of sturmint")

# Check that all required targets are available.
foreach(build ${DRB_BUILD_TYPES})
	if(NOT TARGET "${sturmint_${build}_TARGET}")
		message(FATAL_ERROR "We could not find a ${build} version of sturmint at this location. \
		Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
		rebuild sturmint with a ${build} version as well.")
	endif()
endforeach()

message(STATUS "Found sturmint config at ${sturmint_CONFIG}")
