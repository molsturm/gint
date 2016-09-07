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

cmake_minimum_required(VERSION 3.0.0)

# The DebugReleaseBuild module provides a nice infrastructure to 
# simultaneously build a Debug and a Release version of the same program
# or library, exactly what we need. Try to load it:
include("${CMAKE_SOURCE_DIR}/cmake/IncludeLinalgwrapCmakeModule.cmake")
include_linalgwrap_cmake_module(DebugReleaseBuild)

# initialise DebugReleaseBuild:
drb_init()

##########################################################################
# Global options

# Path where cmake wants the modules and config packages:
set(PackageModuleLocation share/cmake)

# setup global options (needs to be done before call to project())
option(GINT_ENABLE_TESTS "Build gint tests" ON)

##########################################################################
# Project setup

# enable the new version feature introduced in cmake 3.0.0
cmake_policy(SET CMP0048 NEW)
project(gint VERSION 0.0.0)

# Global include directories for this project:
include_directories("${gint_SOURCE_DIR}/src")
include_directories("${gint_BINARY_DIR}/src")

# enable testing of this project
enable_testing()

##########################################################################
# Setup hard and optional dependencies and find components

# Determine and setup compiler flags
drb_setup_compiler_flags()

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

##########################################################################
# Pass on to subdirectories:

# The src directory does most of the work: Here the library is build and 
# installation is configured.
add_subdirectory(src)

# Add subdirectories for the test.
if(GINT_ENABLE_TESTS) 
	add_subdirectory(tests)
endif()

##########################################################################
# Finally install package:
include(cmake/install_package.cmake)
