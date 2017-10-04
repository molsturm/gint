#!/usr/bin/env python3
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
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab

import numpy as np
import gint.gaussian
import gint.sturmian.atomic
from ._available_basis_types import available_basis_types


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


def split_basis_type(basis_type):
    """
    Split the basis type into the type of the
    basis class and the backend specification,
    i.e. returns a tuple (type, string) or
    (type, None) if no backend is specified.

    For example "gaussian" would return (gint.gaussian.Basis, None)
    and "sturmian/atomic/cs_dummy" would return
    (sturmian.atomic.atomic.Basis, "cs_dummy")
    """
    mapping = {
        "gaussian": gint.gaussian.Basis,
        "sturmian/atomic": gint.sturmian.atomic.Basis,
    }

    for key in mapping:
        if basis_type.startswith(key):
            backend = None
            if basis_type[len(key) + 1:]:
                if basis_type not in available_basis_types:
                    raise ValueError("Basis type not available: " + basis_type)
                # Remove the key as well as the trailling "/"
                backend = basis_type[len(key) + 1:]
            return (mapping[key], backend)
    raise ValueError("Unknown basis function type: " + basis_type)
