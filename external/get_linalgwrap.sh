#!/bin/sh
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
## This file is automatically generated from a file in the repository "gscf".
## Edit the original and call the script "update_from_sister_repos.sh" instead.
#
# Settings
#
# The git repo to checkout
FROM="https://github.com/linalgwrap/linalgwrap.git"

# Branch to checkout (empty for default)
BRANCH=""

# Folder to check it out to
WHAT="linalgwrap"

# Interval: How often to update (in seconds)
INTERVAL="3600" # == 1 hour

# File to use in order to test a successful checkout
CHECKFILE="CMakeLists.txt"

if [ ! -f "$PWD/get.lib.sh" ]; then
	echo "PWD needs to be the location of the get.lib.sh file."
	exit 1
fi

. "$PWD/get.lib.sh" || exit 1
update_repo "$FROM" "$WHAT" "$CHECKFILE" "$INTERVAL" "$BRANCH"
exit $?
