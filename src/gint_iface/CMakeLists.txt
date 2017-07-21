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

set(SRC
        construct_structure.cc
        construct_gaussian_basis.cc
)

#
# Parse build type
#
set(GINT_IFACE_BUILD_TYPE "Release" CACHE STRING
"The build type to use for the gint python bindings. \
Valid are Debug and Release.")

set(module_BUILD_TYPE "DEBUG")
if (GINT_IFACE_BUILD_TYPE STREQUAL "Release")
	set(module_BUILD_TYPE "RELEASE")
elseif(GINT_IFACE_BUILD_TYPE  STREQUAL "Debug")
	set(module_BUILD_TYPE "DEBUG")
else()
	message(FATAL_ERROR "The value of GINT_IFACE_BUILD_TYPE needs to be either \"Debug\" or \"Release\".")
endif()

#
# Find swig and python
#
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(Python_ADDITIONAL_VERSIONS 3)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

#
# Download the numpy.i module we require from github
#
file(DOWNLOAD https://raw.githubusercontent.com/numpy/numpy/v1.13.1/tools/swig/numpy.i
	${CMAKE_CURRENT_BINARY_DIR}/numpy.i
	EXPECTED_HASH
	SHA256=bc658c3e9f4ecc5343041339c0697314398fd7cd0e442aa2d424b71086e305ec
)

#
# Setup cpp flags and includes for swig
# (Note, that swig makes heavy use of old-style casts)
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${module_BUILD_TYPE}}")
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-missing-field-initializers)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-pedantic)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-old-style-cast)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-sign-compare)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-maybe-uninitialized)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-deprecated)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-error)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#
# gint swig module
#
set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(gint_iface.i PROPERTIES CPLUSPLUS ON)
swig_add_module(gint_iface python gint_iface.i ${SRC})

swig_link_libraries(gint_iface
	${PYTHON_LIBRARIES}
	${gint_${module_BUILD_TYPE}_TARGET}
)
target_compile_definitions(${SWIG_MODULE_gint_iface_REAL_NAME} PRIVATE
	${module_BUILD_TYPE}=${module_BUILD_TYPE}
)

