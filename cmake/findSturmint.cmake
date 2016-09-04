# Finds and sets up linalgwrap under the target names stored in
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
string(TOUPPER "${PROJECT_NAME}" PROJECT_UPPER)
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

#TODO test and remove (same below in #old#)
message(WARNING "This part of findLinalgwrap has never been tested.")

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
