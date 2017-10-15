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

####################
#-- Empty it all --#
####################
set(GINT_DEPENDENCIES "")
set(GINT_DEPENDENCIES_DEBUG "")
set(GINT_DEPENDENCIES_RELEASE "")
set(GINT_DEPENDENCIES_TEST "")

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
include_krims_cmake_module(FindPackageAutocheckoutFallback)
find_package_autocheckout_fallback(krims 0.2.1)
foreach (build ${DRB_BUILD_TYPES})
	set(GINT_DEPENDENCIES_${build} ${GINT_DEPENDENCIES_${build}} ${krims_${build}_TARGET})
endforeach()

###############
#-- lazyten --#
###############
find_package_autocheckout_fallback(lazyten 0.4.1)
foreach (build ${DRB_BUILD_TYPES})
	set(GINT_DEPENDENCIES_${build} ${GINT_DEPENDENCIES_${build}} ${lazyten_${build}_TARGET})
endforeach()
