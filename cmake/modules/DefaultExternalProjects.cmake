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

# Set cmake policy for ExternalProject module and include it
if (CMAKE_VERSION VERSION_GREATER 3.1)
	cmake_policy(SET CMP0054 NEW)
endif()
include(ExternalProject)
include(CMakeParseArguments)

function(setup_autotools_project NAME)
	# Setup an external project which is configured and built
	# using autotools.
	#
	# Syntax:
	#   setup_autotools_project(project
	#          GIT_REPOSITORY "https://gitserver/repo.git"
	#          CONFIGURE_OPTS "--with-blubber"
	#   )
	#
	# Basically the name is the first argument, followed by the usual options
	# to ExternalProject_Add and then as a last argument CONFIGURE_OPTS.
	#
	# Note: This order is important!
	#

	set(options )
	set(oneValueArgs )
	set(multiValueArgs CONFIGURE_OPTS)
	cmake_parse_arguments(SAP "${options}" "${oneValueArgs}" "${multiValueArgs}"  ${ARGN})

	# Disable automatic updates in case of cmake > 3.2
	if (CMAKE_VERSION VERSION_GREATER 3.2)
		set(UPDATE_OPTIONS UPDATE_DISCONNECTED ON)
	endif()

	set(configure_wrap "${PROJECT_BINARY_DIR}/scripts/autotools_configure.sh")
	file(WRITE ${configure_wrap} 
"#!/bin/sh

# Script which calls configure in the current directory
# after this script has been generated using the appropriate
# toos from the autotools suite.

# ------------------------------------------------

are_autotools_available() {
	which autoconf > /dev/null 2>&1 || return 1
	return 0
}

# ------------------------------------------------

if ! are_autotools_available; then
	echo \"Some required programs of the autotools suite could not be found.\" >&2
	exit 1
fi

SOURCE_DIR=\"$1\"
shift

if [ ! -d \"$SOURCE_DIR\" ]; then
	echo \"Could not find source directory $SOURCE_DIR\" >&2
	exit 1
fi
if [ ! -x \"$SOURCE_DIR/autogen.sh\" ]; then
	echo \"Could not find autogen.sh script inside source directory $SOURCE_DIR.\" >&2
	echo \"This script is needed to generate the configure.sh script\" >&2
	exit 1
fi

if [ ! -f \"$SOURCE_DIR/configure.ac\" ]; then
	echo \"Could not find configure.ac inside source directory $SOURCE_DIR.\" >&2
	exit 1
fi


# The configure script needs an update or does not exist:
if [ ! -f \"$SOURCE_DIR/configure\" -o \"$SOURCE_DIR/configure.ac\" -nt \"$SOURCE_DIR/configure\" ]; then
	OLDDIR=`pwd`
	cd \"$SOURCE_DIR\"
	./autogen.sh
	cd \"$OLDDIR\"
fi

# Configure needs to be re-run
if [ ! -f config.status -o  \"$SOURCE_DIR/configure\" -nt config.status ]; then
	$SOURCE_DIR/configure \"$@\"
	exit $?
fi

exit 0
"	)

	# TODO The problem with this behaviour is that it can lead to double-parallelisation:
	#      Both the building of the external project and of the outer cmake build try to access
	#      all cpus.
	option(${NAME}_BUILD_PARALLEL "Build external integral libraries in parallel." ON)

	# Determine the processor count for parallel build
	include(ProcessorCount)
	ProcessorCount(NCPU)
	if (NOT NCPU EQUAL 0 AND ${NAME}_BUILD_PARALLEL)
		set(MAKE_ARGS -j ${NCPU})
	endif()

	ExternalProject_Add(${NAME}
		PREFIX "${PROJECT_BINARY_DIR}/external/${NAME}"
		${UPDATE_OPTIONS}
		#
		CONFIGURE_COMMAND /bin/sh ${configure_wrap} <SOURCE_DIR> --prefix=<INSTALL_DIR> ${SAP_CONFIGURE_OPTS}
		BUILD_COMMAND make ${MAKE_ARGS}
		INSTALL_COMMAND make ${MAKE_ARGS} install
		#
		${SAP_UNPARSED_ARGUMENTS}
	)
endfunction()

