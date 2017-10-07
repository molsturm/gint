#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the sturmint authors
##
## This file is part of sturmint.
##
## sturmint is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## sturmint is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with sturmint. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

import numpy as np


def cartesian_to_spherical(x, y, z):
    """
    Return the tuple r, theta, phi, which gives the points x, y, z
    in equivalent radial corrdinates.
    """
    xy_sq = x**2 + y**2
    r = np.sqrt(xy_sq + z**2)
    theta = np.arctan2(np.sqrt(xy_sq), z)
    phi = np.arctan2(y, x)

    # Note: Using the trigonometric identities one can show that
    #
    #       np.arctan2(np.sqrt(xy_sq), z) == np.arccos( z / r)
    #
    #       The former expression is better since it avoids the costly
    #       division and allows better vectorisation (since only one type
    #       of function used.
    return r, theta, phi
