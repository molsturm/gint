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

cmake_minimum_required(VERSION 3.0.0)
cmake_policy(VERSION 3.0.0)

# Load the DebugReleaseBuild module for a simultaneous Debug and Release
# and release build and initialise it.
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/IncludeKrimsCmakeModule.cmake")
include_krims_cmake_module(DebugReleaseBuild)
drb_init()

##########################################################################
# Project setup

project(gint VERSION 0.0.0)
message("\n# Configuring ${PROJECT_NAME}\n#")

# Global include directories for this project:
include_directories("${gint_SOURCE_DIR}/src")
include_directories("${gint_BINARY_DIR}/src")

# Extra module path for this project:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

# enable testing of this project
enable_testing()

##########################################################################
# Global options

include_krims_cmake_module(GlobalOption)
global_option(ENABLE_DOCUMENTATION "Build documentation"                       OFF )
global_option(ENABLE_EXAMPLES      "Build example executables"                 ON  )
global_option(ENABLE_TESTS         "Build unit test executables"               ON  )
global_option(INTERFACE_PYTHON     "Build and install python interface module" ON  )

##########################################################################
# Setup hard and optional dependencies and find components

# Determine and setup compiler flags, enforcing C++11
drb_setup_compiler_flags(11)

# Export compile commands if cmake supports it.
if (CMAKE_VERSION VERSION_GREATER 3.5.0)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# If a user wants a python module, setup the interface helper class
# This might alter the compiler flags in order for the compiled
# code to be compatible with the interfacing, hence this call is
# needed right here.
if (GINT_INTERFACE_PYTHON)
	include_krims_cmake_module(CppInterface)
endif()

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

# Setup download and installation of static data
include(cmake/static_data.cmake)

##########################################################################
# Pass on to subdirectories:

# The src directory does most of the work: Here the library is build and
# installation is configured.
add_subdirectory(src)
set(GINT_SUBDIRS src)

# Add subdirectories for the testing.
if(GINT_ENABLE_TESTS)
	add_subdirectory(tests)
	set(GINT_SUBDIRS ${GINT_SUBDIRS} tests)
endif()

# Add subdirectories for the examples.
if(GINT_ENABLE_EXAMPLES)
	add_subdirectory(examples)
	set(GINT_SUBDIRS ${GINT_SUBDIRS} examples)
endif()


if (GINT_ENABLE_DOCUMENTATION)
	add_subdirectory(doc)
endif()

##########################################################################
# Setup extra targets for clang tools
include_krims_cmake_module(SetupClangTargets)
add_available_clang_targets_for(${PROJECT_NAME}
	DIRECTORIES ${GINT_SUBDIRS})

# Finally install package:
include(cmake/install_package.cmake)

# and print a closing message
message("#\n# Configuring ${PROJECT_NAME} done\n")
