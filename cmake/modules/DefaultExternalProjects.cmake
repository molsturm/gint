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
	#          CHECKPOINT_FROM "fromfile"
	#          CHECKPOINT_TO "fromfile"
	#          CONFIGURE_OPTS "--with-blubber"
	#   )
	#
	# Basically the name is the first argument, followed by the usual options
	# to ExternalProject_Add and then as a last argument CONFIGURE_OPTS.
	#
	# The other options allow to define two checkpoint files. Make install
	# will only be called if CHECKPOINT_FROM is newer than CHECKPOINT_TO
	# CHECKPOINT_TO is also added as a byproduct to the ExternalProject_Add
	# call if no byproducts are specified.
	#
	# Note: This order is important!
	#

	set(options)
	set(oneValueArgs CHECKPOINT_TO CHECKPOINT_FROM)
	set(multiValueArgs CONFIGURE_OPTS BUILD_BYPRODUCTS)
	cmake_parse_arguments(SAP "${options}" "${oneValueArgs}" "${multiValueArgs}"  ${ARGN})

	# Set byproducts if not already specified
	if ("${SAP_BUILD_BYPRODUCTS}" STREQUAL "")
		set(SAP_BUILD_BYPRODUCTS "${SAP_CHECKPOINT_TO}")
	endif()

	if (NOT "${SAP_BUILD_BYPRODUCTS}" STREQUAL "")
		if (CMAKE_VERSION VERSION_GREATER 3.2)
			# The BUILD_BYPRODUCTS flag is only understood from 3.2 on
			set(BUILD_BYPRODUCT_OPTIONS
				"BUILD_BYPRODUCTS" ${SAP_BUILD_BYPRODUCTS})
		endif()
	elseif(CMAKE_GENERATOR STREQUAL "Ninja")
		message(FATAL_ERROR "BUILD_BYPRODUCTS is required if Ninja is used as a generator! \
Either specify CHECKPOINT_TO or BUILD_BYPRODUCTS in the call")
	endif()

	# Set to the default special values, which mean that this feature is ignored
	if ("${SAP_CHECKPOINT_TO}" STREQUAL "")
		set(SAP_CHECKPOINT_TO "<none>")
	endif()
	if ("${SAP_CHECKPOINT_FROM}" STREQUAL "")
		set(SAP_CHECKPOINT_FROM "<none>")
	endif()

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
	echo \"Running autogen.sh to generate configure script.\"
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

	set(install_wrap "${PROJECT_BINARY_DIR}/scripts/autotools_install.sh")
	file(WRITE ${install_wrap}
"#!/bin/sh

# Script which calls make install in the current directory
# but only if a certain checkpoint file is not older than the
# last time make install was run.

# ------------------------------------------------

CHECKPOINT_FROM=\"$1\"
CHECKPOINT_TO=\"$2\"
shift
shift

FORCE_INSTALL=0
if [ \"$CHECKPOINT_FROM\" = \"<none>\" ] && [ \"$CHECKPOINT_TO\" = \"<none>\" ]; then
	FORCE_INSTALL=1
fi

if [ ! -f \"$CHECKPOINT_FROM\" ]; then
	echo \"Could not find checkpoint from file $CHECKPOINT_FROM\" >&2
	exit 1
fi

if [ \"x$FORCE_INSTALL\" = \"x1\" ] || [ ! -f \"$CHECKPOINT_TO\" ] || [ \"$CHECKPOINT_FROM\" -nt \"$CHECKPOINT_TO\" ]; then
	make $@ install
	exit $?
fi

exit 0
"	)


	# TODO The problem with this behaviour is that it can lead to double-parallelisation:
	#      Both the building of the external project and of the outer cmake build try to access
	#      all cpus.
	include(ProcessorCount)
	ProcessorCount(NCPU)
	set(DEP_BUILD_${NAME}_NJOBS ${NCPU} CACHE STRING
		"Build external project ${NAME} in parallel with this number of jobs (1 or 0 disables parallel build)")

	if ("${DEP_BUILD_${NAME}_NJOBS}" GREATER 1)
		set(MAKE_ARGS -j ${DEP_BUILD_${NAME}_NJOBS})
	endif()

	ExternalProject_Add(${NAME}
		PREFIX "${PROJECT_BINARY_DIR}/external/${NAME}"
		${UPDATE_OPTIONS}
		${BUILD_BYPRODUCT_OPTIONS}
		#
		CONFIGURE_COMMAND /bin/sh ${configure_wrap} <SOURCE_DIR> --prefix=<INSTALL_DIR> ${SAP_CONFIGURE_OPTS}
		BUILD_COMMAND make ${MAKE_ARGS}
		INSTALL_COMMAND /bin/sh ${install_wrap} ${SAP_CHECKPOINT_FROM} ${SAP_CHECKPOINT_TO} ${MAKE_ARGS} install
		#
		${SAP_UNPARSED_ARGUMENTS}
	)
endfunction()

