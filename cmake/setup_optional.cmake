# Setup optional dependencies and features
# alters these things
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

##########################
#--  Static integrals  --#
##########################
option(GINT_ENABLE_STATIC_INTEGRALS "Enable a basis types which consist entirely of pre-computed integral data. \
(Increases binary size, but useful for testing)" OFF)

if(GINT_ENABLE_STATIC_INTEGRALS)
	message(STATUS "Enabled pre-computed static integrals")
	LIST(APPEND GINT_DEFINITIONS "GINT_HAVE_STATIC_INTEGRALS")
endif()
