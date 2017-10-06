#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the molsturm authors
##
## This file is part of molsturm.
##
## molsturm is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## molsturm is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with molsturm. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

from .._available_basis_types import available_basis_types
from . import _gaussian_shells as _expr
from .. import _iface
from ..Structure import Structure
import collections
import numpy as np


class Shell(collections.namedtuple("Shell", ["l", "pure", "contraction_coefficients",
                                             "contraction_exponents", "centre"])):
    @property
    def n_bas(self):
        return 2 * self.l + 1 if self.pure else (self.l + 1) * (self.l + 2) // 2

    def __len__(self):
        return self.n_bas()


"""The list of available gaussian backends"""
available_backends = [t[9:] for t in available_basis_types if t.startswith("gaussian/")]


class Basis:
    def __init__(self, structure, basis_set_name, backend="auto"):
        """
        structure    Either a gint.Structure object or a tuple (atoms, coords)
                     or just a single atom.
                     which gives the list of atoms and the list of their coordinates
        basis_set_name   String describing the gaussian basis set to use
        backend      Specify the precise backend to use.
                     Should be available. Otherwise an error is produced.
        """
        if len(available_backends) == 0:
            raise RuntimeError("No gaussian backend is available in gint. The list of "
                               "available basis types is: " +
                               ",".join(available_basis_types))

        if not isinstance(structure, Structure):
            if isinstance(structure, (str, int)):
                structure = Structure(atoms=structure)
            else:
                structure = Structure(*structure)

        basis_raw = _iface.construct_gaussian_basis(structure.atom_numbers,
                                                    structure.coords, basis_set_name)
        shells_raw = (basis_raw.shell(i) for i in range(basis_raw.n_shells()))
        self.basis_set_name = basis_set_name

        # Setup the shells.
        # Note the explicit copy! This is needed, since the basis_raw.shell
        # function merely returns a view into the data, which is held by
        # the c++ object behing basis_raw.
        self.shells = [Shell(sh.l, sh.pure, sh.coefficients().copy(),
                             sh.exponents().copy(), sh.origin().copy())
                       for sh in shells_raw]

        if backend == "auto" or backend is None:
            # List the priority of the backends
            __backend_priority = ["libint"]

            for b in __backend_priority + available_backends:
                if b in available_backends:
                    self.backend = b
                    break
        else:
            if backend not in available_backends:
                raise ValueError("The gaussian inegral backend " + backend + " is not "
                                 "available. The following basis types are implemented:"
                                 " " + ",".join(available_basis_types))
            self.backend = backend

        if self.backend == "libint":
            # Libint as we configure it uses standard ordering,
            # i.e. the ordering of the Common Component Architecture.
            # (xxx, xxy, xxz, xyy, xyz, xzz, yyy, ...)
            self.cartesian_ordering = "standard"

    @property
    def basis_type(self):
        return "gaussian/" + self.backend

    def __len__(self):
        return sum(s.n_bas for s in self.shells)

    @property
    def size(self):
        """Return the number of basis functions in this basis"""
        return self.__len__()

    @property
    def n_bas(self):
        """Return the number of basis functions in this basis"""
        return self.__len__()

    @property
    def n_shells(self):
        """Return the number of shells in this basis"""
        return len(self.shells)

    @property
    def has_real_harmonics(self):
        """
        Does this basis have real functions to describe the angular part. This impies that
        the repulsion integrals satisfy the extra symmetry (ab|cd) = (ba|cd) in shell
        pair notation, which would not be true otherwise
        """
        return True

    def evaluate_at(self, x, y, z, mask=None):
        if mask:
            raise NotImplementedError("Basis mask is not yet implemented.")

        evaluated = []
        for sh in self.shells:
            if sh.pure:
                res = _expr.pure_shell(x, y, z, sh.l, sh.centre,
                                       sh.contraction_coefficients,
                                       sh.contraction_exponents)
            else:
                res = _expr.cartesian_shell(x, y, z, sh.l, sh.centre,
                                            sh.contraction_coefficients,
                                            sh.contraction_exponents,
                                            self.cartesian_ordering)
            evaluated.append(res)
        return np.concatenate(evaluated)
