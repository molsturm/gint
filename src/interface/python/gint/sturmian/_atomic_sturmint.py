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

from ..BasisBase import BasisBase, available_basis_types
import sturmint.atomic.cs

"""The list of available sturmian backends"""
available_backends = [t[16:] for t in available_basis_types
                      if t.startswith("sturmian/atomic")]


class Basis(sturmint.atomic.cs.Basis, BasisBase):
    def __init__(self, structure, k_exp, n_max, l_max=None, m_max=None, backend="auto"):
        if backend == "auto" or backend is None:
            # List the priority of the backends
            __backend_priority = ["cs_reference_pc", "cs_dummy", "cs_static14"]

            for b in __backend_priority + available_backends:
                if b in available_backends:
                    self.backend = b
                    break
        else:
            if backend not in available_backends:
                raise ValueError("The Coulomb-Sturmian integral backend " + backend +
                                 " is not available. The following basis types "
                                 "are implemented: " +
                                 ",".join(available_basis_types))
            self.backend = backend

        # TODO If the backend does not use nlm order, change that here.
        super().__init__(k_exp, n_max, l_max, m_max, order="nlm")

    @property
    def has_real_harmonics(self):
        """
        Does this basis have real functions to describe the angular part. This impies that
        the repulsion integrals satisfy the extra symmetry (ab|cd) = (ba|cd) in shell
        pair notation, which would not be true otherwise
        """
        # Currently all backends use complex spherical harmonics, but for m_max == 0 it
        # does not matter, since the real and complex harmonics are the same, so we
        # prefer real in that case.
        return all(nlm[2] <= 0 for nlm in self.functions)

    @property
    def size(self):
        """Return the number of basis functions in this basis"""
        return self.__len__()

    @property
    def basis_type(self):
        return "sturmian/atomic/" + self.backend
