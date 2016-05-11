# Finds and sets up the sturmint targets
#
#   Upstream::sturmint.g      (Debug target)
#   Upstream::sturmint        (Release target)
#
# Such that just linking against them as a dependency does everything
# automaatically.
#
# If sturmint was not compiled as Debug and Release it may happen that
# one of the above targets is not available.
#
# For convenience and future portability the variables 
#      sturmint_DEBUG_TARGET      (Name of Debug target)
#      sturmint_RELEASE_TARGET    (Name of Release target)
# are available, but note that their presence does not indicate that
# the targets are actually available.
#

# We need sturmint version 0.0.0
find_package(sturmint 0.0.0 REQUIRED CONFIG 
	PATHS 
	${CMAKE_SOURCE_DIR}/../sturmint/build
)
# Now we found the library. Most of the times that's it and we are done.
# But if we got the sturmint from a build directory, then it is very
# likely that the header includes cannot be found like this.

# Find dir containing sturmint config:
get_filename_component(sturmint_config_dir "${sturmint_CONFIG}"  DIRECTORY)

set(sturmint_DEBUG_TARGET   "Upstream::sturmint.g" 
	CACHE INTERNAL "Target name of debug version of sturmint")
set(sturmint_RELEASE_TARGET "Upstream::sturmint"
	CACHE INTERNAL "Target name of release version of sturmint")

# So we try to correct this for all existing targets:
foreach(target ${sturmint_DEBUG_TARGET} ${sturmint_RELEASE_TARGET})
	if (NOT TARGET ${target})
		continue()
	endif()

	# Check that the includes are present:
	get_target_property(STURMINT_INTERFACE_INCLUDES 
		${target} INTERFACE_INCLUDE_DIRECTORIES)

	# If yes continue
	if(NOT "${STURMINT_INTERFACE_INCLUDES}" 
			STREQUAL "STURMINT_INTERFACE_INCLUDES-NOTFOUND")
		continue()
	endif()

	set(STURMINT_INTERFACE_INCLUDES "")

	# Try to find the include dirctory:
	find_path(sturmint_INCLUDE_DIR "sturmint/version.hh"
		HINTS
		# If we found a build directory, then this is the 
		# path to the include directory
		${sturmint_config_dir}/../src/
		PATHS
		$ENV{sturmint_INCLUDE_DIR}
		${CMAKE_SOURCE_DIR}/../sturmint/src
		DOC "sturmint header include directory"
	)

	# Check that the include directory was found
	if ("${sturmint_INCLUDE_DIR}" STREQUAL "sturmint_INCLUDE_DIR-NOTFOUND")
		message(FATAL_ERROR "Could not find sturmint include directory. 
Please provide a hint using the environment variable sturmint_INCLUDE_DIR")
	endif()

	# Append to interface includes:
	set(STURMINT_INTERFACE_INCLUDES ${STURMINT_INTERFACE_INCLUDES} ${sturmint_INCLUDE_DIR})

	# Check that the sturmint/version_defs.hh file can be found in this include
	# directory.
	if(NOT EXISTS "${sturmint_INCLUDE_DIR}/sturmint/version_defs.hh")
		if(EXISTS "${sturmint_config_dir}/src/sturmint/version_defs.hh")
			set(STURMINT_INTERFACE_INCLUDES ${STURMINT_INTERFACE_INCLUDES} 
				"${sturmint_config_dir}/src"
			)
		else()
			message(FATAL_ERROR "Could not find sturmint version_defs.hh file")
		endif()
	endif()

	# Set the interface includes:
	set_target_properties(${target}
		PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${STURMINT_INTERFACE_INCLUDES}"
	)
endforeach()

message(STATUS "Found sturmint config at ${sturmint_CONFIG}")
