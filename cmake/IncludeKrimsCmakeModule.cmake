## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the gint authors
##
## This file is part of gint.
##
## gint is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gint is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with gint. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

## DO NOT EDIT
## This file is automatically generated from a file in the repository "linalgwrap".
## Edit the original and call the script "update_from_sister_repos.sh" instead.

# Try to find a cmake module which gets shipped with krims.
#
# In case the module is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, krims is automatically checked out and the module is loaded from there.
# Otherwise a fatal error is produced.
#

option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)

if (NOT DEFINED EXTERNAL_DIR)
	if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/../external/get_krims.sh")
		set(EXTERNAL_DIR "${CMAKE_CURRENT_LIST_DIR}/../external")
	else()
		message(FATAL_ERROR "Could not determine the location of the external \
projects. Try to specify it manually using the variable EXTERNAL_DIR.")
	endif()
endif()

function(include_krims_cmake_module MODULE)
	# Sanity check
	if(NOT EXISTS "${EXTERNAL_DIR}/get_krims.sh")
		message(FATAL_ERROR "get_krims.sh does not exist in EXTERNAL_DIR \
== ${EXTERNAL_DIR}.")
	endif()

	# First try to load it plainly as a module:
	include(${MODULE} OPTIONAL RESULT_VARIABLE RES)
	if (NOT "${RES}" STREQUAL "NOTFOUND")
		return()
	endif()

	# We could not "just" find it. Try the krims_DIR hint:
	include("$ENV{krims_DIR}/share/cmake/modules/${MODULE}.cmake"
		OPTIONAL RESULT_VARIABLE RES)
	if (NOT "${RES}" STREQUAL "NOTFOUND")
		return()
	endif()

	# Try autocheckout
	if (AUTOCHECKOUT_MISSING_REPOS)
		execute_process(
			COMMAND "sh" "get_krims.sh"
			WORKING_DIRECTORY "${EXTERNAL_DIR}"
			RESULT_VARIABLE RES
			)
		if (NOT RES EQUAL 0)
			message(FATAL_ERROR "Getting krims from git failed with error: ${RES}")
		endif()

		include("${EXTERNAL_DIR}/krims/cmake/modules/${MODULE}.cmake")
		return()
	endif()

	message(FATAL_ERROR "Could not find the ${MODULE} module.
Either provide the installation prefix of the krims library in the environment \
variable krims_DIR or enable autocheckout via '-DAUTOCHECKOUT_MISSING_REPOS=ON'.")
endfunction(include_krims_cmake_module)
