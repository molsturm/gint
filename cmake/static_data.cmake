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

# This module deals with downloading and installing static data
# which is needed to drive some of the sturmint integrals.

include_krims_cmake_module(DataFiles)

# The name of the target to download the data
set(DATA_DOWN_TARGET "${PROJECT_NAME}-download-data")

# Download the data tarball and unpack it into the data directory
set(TARHASH 40e1d61e111adec23ba4d05d260a79c3e66fa1635eff6980e092ccaf02b83133)
data_download_target(${DATA_DOWN_TARGET}
	"https://get.molsturm.org/data_tar/${PROJECT_NAME}/${TARHASH}.tar.gz"
	SHA256=${TARHASH}
)
unset(TARHASH)

# Install the data to the ${PROJECT_NAME}_DATA_INSTALL_DIR
install_data()
