#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
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

from .BasisBase import available_basis_types
import gint.gaussian
import gint.sturmian.atomic


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


class Basis:
    """
    Common class to construct arbitrary gint bases.
    """

    @classmethod
    def construct(cls, basis_type, structure, **kwargs):
        """
        Construct any of the available gint basis objects using
        the basis type string and the arguments needed to construct
        the basis object.

        basis_type:   String identifying the basis type
        structure:    The gint.Structure object of the structure to model

        Example:
            struct = gint.Structure("Be")
            Basis.construct("gaussian/libint", struct, basis_set_name="sto-3g")
        """
        Basis, backend = split_basis_type(basis_type)
        try:
            return Basis(structure, **kwargs, backend=backend)
        except (TypeError, ValueError) as e:
            raise ValueError("Invalid argument for basis construction: " + str(e))
