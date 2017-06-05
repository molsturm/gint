# Setup optional dependencies and features
# alters these things
#
# 	GINT_DEPENDENCIES		everyone needs these libraries
# 	GINT_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	GINT_DEPENDENCIES_RELEASE	release mode needs these extras
# 	GINT_DEPENDENCIES_TEST		tests need these extra libraries
#
#       GINT_DEFINITIONS		definitions for all compilation
#       GINT_DEFINITIONS_DEBUG		definitions for debug mode
#       GINT_DEFINITIONS_RELEASE	definitions for release mode
#

####################
#-- C++ standard --#
####################
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 14)
	message(STATUS "Detected C++14 support: Setting GINT_HAVE_CXX14")
	LIST(APPEND GINT_DEFINITIONS "GINT_HAVE_CXX14")
endif()
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 17)
	message(STATUS "Detected C++17 support: Setting GINT_HAVE_CXX17")
	LIST(APPEND GINT_DEFINITIONS "GINT_HAVE_CXX17")
endif()

################
#--  Libint2 --#
################
option(GINT_ENABLE_LIBINT "Enable the libint library to compute Gaussian integrals." OFF)
option(GINT_LIBINT_USE_SYSTEM "Enable the use of a system-provided libint library" OFF)
set(GINT_LIBINT_MAX_AM 6 CACHE STRING
"Maximal angular momentum libint can perform integrals over. \
Choose a smaller value to get a faster build."
)   # Note: ORCA uses a value of 7 (up to K) in the setting above.

if (GINT_ENABLE_LIBINT)
	# Check options:
	if (GINT_LIBINT_MAX_AM LESS 4)
		message(FATAL_ERROR "GINT_LIBINT_MAX_AM needs to be at least 4")
	endif()

	# Forward parameters to included module
	set(LIBINT_VERSION 2.3.1) # We need at least this version
	set(LIBINT_SEARCH_SYSTEM ${GINT_LIBINT_USE_SYSTEM})
	set(LIBINT_MAX_AM ${GINT_LIBINT_MAX_AM})

	include(cmake/findLibint.cmake)
	unset(LIBINT_VERSION)
	unset(LIBINT_SEARCH_SYSTEM)

	LIST(APPEND GINT_DEFINITIONS "GINT_HAVE_LIBINT")
	set(GINT_DEPENDENCIES ${LIBINT_TARGET})
endif()

##########################
#--  Static integrals  --#
##########################
option(GINT_ENABLE_STATIC_INTEGRALS "Enable a basis types which consist entirely of pre-computed integral data. \
(Increases binary size, but useful for testing)" OFF)

if(GINT_ENABLE_STATIC_INTEGRALS)
	message(STATUS "Enabled pre-computed static integrals")
	LIST(APPEND GINT_DEFINITIONS "GINT_HAVE_STATIC_INTEGRALS")
endif()
