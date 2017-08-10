#!/usr/bin/env python3
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab
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

from .. import element
from . import _gaussian_shells as _expr
from .. import _iface
from ..interface import Structure
import collections
import numpy as np

class Shell(collections.namedtuple("Shell", [ "l", "pure", "contraction_coefficients",
                                              "contraction_exponents", "centre" ])):
  @property
  def n_bas(self):
    return 2 * self.l + 1 if self.pure else (self.l + 1)*(self.l + 2) // 2
  def __len__(self):
    return self.n_bas()


class Basis:
  def __init__(self, structure, basis_set, basis_type=None,
               cartesian_ordering="standard"):
    """
    structure    Either a gint.Structure object or a tuple (atoms, coords)
                 or just a single atom.
                 which gives the list of atoms and the list of their coordinates
    basis_set    String describing the gaussian basis set to use
    ordering     The ordering of the basis functions to use if cartesian
                 angular functions are used.
                 For pure shells (i.e. with cartesian spherical harmonics)
                 the ordering is from m=-2 (The terms with more sines) to
                 m=+2 (The terms with more coses)

                 The currently implemented orderings are
                   - standard:  The ordering of the Common Component Architecture.
                                This is just (xxx, xxy, xxz, xyy, xyz, xzz, yyy, ...)
    basis_type   Specify a basis type as well (optional)
    """
    if not isinstance(structure, Structure):
      if isinstance(structure, (str,int)):
        structure = Structure(atoms=structure)
      else:
        structure = Structure(*structure)

    basis_raw = _iface.construct_gaussian_basis(structure.atom_numbers, structure.coords,
                                                basis_set)
    shells_raw = ( basis_raw.shell(i) for i in range(basis_raw.n_shells()) )

    # Setup the shells.
    # Note the explicit copy! This is needed, since the basis_raw.shell
    # function merely returns a view into the data, which is held by
    # the c++ object behing basis_raw.
    self.shells = [ Shell(sh.l, sh.pure, sh.coefficients().copy(),
                          sh.exponents().copy(), sh.origin().copy())
                    for sh in shells_raw ]

    self.cartesian_ordering = cartesian_ordering
    if cartesian_ordering != "standard":
      raise ValueError("Currently only standard cartesian_ordering is implemented.")

    self.basis_type = basis_type

  def evaluate_at(self, x, y, z, mask=None):
    if mask:
      raise NotImplementedError("Basis mask is not yet implemented.")

    evaluated=[]
    for sh in self.shells:
      if sh.pure:
        res = _expr.pure_shell(x, y, z, sh.l, sh.centre,
                               sh.contraction_coefficients, sh.contraction_exponents)
      else:
        res = _expr.cartesian_shell(x, y, z, sh.l, sh.centre,
                                    sh.contraction_coefficients,
                                    sh.contraction_exponents,
                                    self.cartesian_ordering)
      evaluated.append(res)
    return np.concatenate(evaluated)

