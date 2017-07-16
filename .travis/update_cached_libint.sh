#!/bin/sh
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

# The directory where the libint cached files should live
CACHEDIR="$HOME/cache_libint"

# Place where the libint files are installed by ExternalProject
LIBINT_INSTALL_DIR="${TRAVIS_BUILD_DIR}/build/external/libint"

# -------------------------------------------------------------

# Replace cache by our build.
if [ -d "$LIBINT_INSTALL_DIR" ]; then
	echo "Updating libint cache at $CACHEDIR"

	rm -rf "$CACHEDIR"
	mkdir "$CACHEDIR"

	for subdir in include lib share; do
		cp -a $LIBINT_INSTALL_DIR/$subdir $CACHEDIR
	done
fi

