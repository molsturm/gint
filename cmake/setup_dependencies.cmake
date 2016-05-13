# sets these things
#
# 	GINT_DEPENDENCIES			everyone needs these libraries
# 	GINT_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	GINT_DEPENDENCIES_RELEASE		release mode needs these extras
# 	GINT_DEPENDENCIES_TEST		tests need these extra libraries
#
#       GINT_DEFINITIONS			definitions for all compilation
#       GINT_DEFINITIONS_DEBUG		definitions for debug mode
#       GINT_DEFINITIONS_RELEASE		definitions for release mode
#       

####################
#-- Empty it all --#
####################
set(GINT_DEPENDENCIES "")
set(GINT_DEPENDENCIES_DEBUG "")
set(GINT_DEPENDENCIES_RELEASE "")
set(GINT_DEPENDENCIES_TEST "")
set(GINT_DEFINITIONS "")
set(GINT_DEFINITIONS_DEBUG "")
set(GINT_DEFINITIONS_RELEASE "")

##################
#-- linalgwrap --#
##################
if (TARGET ${linalgwrap_DEBUG_TARGET} OR TARGET ${linalgwrap_RELEASE_TARGET})
	# If the targets are already defined elsewhere, we are done:
	message(STATUS "Using linalgwrap library provided by build environment.")
else()
	include(cmake/findLinalgwrap.cmake)
endif()

# Check that all required targets are available.
foreach(build ${DRB_BUILD_TYPES})
	if(NOT TARGET "${linalgwrap_${build}_TARGET}")
		message(FATAL_ERROR "We could not find a ${build} version of linalwrap at this location. \
Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
rebuild linalgwrap with a ${build} version as well.")
	endif()

	# Add dependencies to appropriate versions of gint
	set(GINT_DEPENDENCIES_${build} ${GINT_DEPENDENCIES_${build}} ${linalgwrap_${build}_TARGET})
endforeach()

################
#-- sturmint --#
################

if (TARGET ${sturmint_DEBUG_TARGET} OR TARGET ${sturmint_RELEASE_TARGET})
	message(STATUS "Using sturmint library provided by build environment.")
else()
	include(cmake/findSturmint.cmake)
endif()

# Check that all required targets are available.
foreach(build ${DRB_BUILD_TYPES})
	if(NOT TARGET "${sturmint_${build}_TARGET}")
		message(FATAL_ERROR "We could not find a ${build} version of sturmint at this location. \
Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
rebuild sturmint with a ${build} version as well.")
	endif()

	# Add dependencies to appropriate versions of gint
	set(GINT_DEPENDENCIES_${build} ${GINT_DEPENDENCIES_${build}} ${sturmint_${build}_TARGET})
endforeach()


##############
#-- catch  --#
##############
if(GINT_ENABLE_TESTS)
	if (TARGET common_catch)
		MESSAGE(STATUS "Using catch provided by build enviroment.")
		set(GINT_DEPENDENCIES_TEST ${GINT_DEPENDENCIES_TEST} common_catch)
	else()
		include(cmake/findCatch.cmake)
		set(GINT_DEPENDENCIES_TEST ${GINT_DEPENDENCIES_TEST} Upstream::catch)
	endif()

endif()

##################
#-- rapidcheck --#
##################
if(GINT_ENABLE_TESTS)
	if (TARGET common_rapidcheck)
		MESSAGE(STATUS "Using rapidcheck provided by build environment.")
		set(GINT_DEPENDENCIES_TEST ${GINT_DEPENDENCIES_TEST} rapidcheck)
	else()
		find_package(rapidcheck REQUIRED CONFIG 
			PATHS 
			${gint_SOURCE_DIR}/../linalgwrap/build
		)

		message(STATUS "Found rapidcheck config at ${rapidcheck_CONFIG}")
		set(GINT_DEPENDENCIES_TEST ${GINT_DEPENDENCIES_TEST} Upstream::rapidcheck)
	endif()
endif()
