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

# sets these things
#
#       GINT_DEPENDENCIES               everyone needs these libraries
#       GINT_DEPENDENCIES_DEBUG         debug mode needs these extras
#       GINT_DEPENDENCIES_RELEASE       release mode needs these extras
#       GINT_DEPENDENCIES_TEST          tests need these extra libraries
#
#       GINT_DEFINITIONS                definitions for all compilation
#       GINT_DEFINITIONS_DEBUG          definitions for debug mode
#       GINT_DEFINITIONS_RELEASE        definitions for release mode
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

###############
#--  Types  --#
###############
# Determine scalar types we want to use here:
include_krims_cmake_module(ScalarTypes)
setup_scalar_types()
set(GINT_DEPENDENCIES ${GINT_DEPENDENCIES} ${SCALAR_TYPES_LIBRARIES})

############################
#-- rapidcheck and catch --#
############################
if(GINT_ENABLE_TESTS)
	# We need to setup rapidcheck and catch for the tests:
	include(cmake/findRapidcheck.cmake)
	set(GINT_DEPENDENCIES_TEST ${GINT_DEPENDENCIES_TEST} ${rapidcheck_TARGET})

	include(cmake/findCatch.cmake)
	set(GINT_DEPENDENCIES_TEST ${GINT_DEPENDENCIES_TEST} ${catch_TARGET})
endif()

#############
#-- krims --#
#############
# Find at least version 0.0.0
set(KRIMS_VERSION 0.0.0)
include(cmake/findKrims.cmake)

foreach (build ${DRB_BUILD_TYPES})
	set(GINT_DEPENDENCIES_${build} ${GINT_DEPENDENCIES_${build}} ${krims_${build}_TARGET})
endforeach()

##################
#-- linalgwrap --#
##################
# Find at least version 0.2.0
set(LINALGWRAP_VERSION 0.2.0)
include(cmake/findLinalgwrap.cmake)

foreach (build ${DRB_BUILD_TYPES})
	set(GINT_DEPENDENCIES_${build} ${GINT_DEPENDENCIES_${build}} ${linalgwrap_${build}_TARGET})
endforeach()
