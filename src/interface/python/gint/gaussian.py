#!/usr/bin/env python3
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
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab

from . import element
import collections
import numpy as np
from . import _iface
from .interface import Structure

Shell = collections.namedtuple("Shell", [ "l", "pure", "contraction_coefficients",
                                          "contraction_exponents", "centre" ])


class Basis:
  def __init__(self, structure, basis_set, ordering="libint"):
    """
    structure    Either a gint.Structure object or a tuple (atoms, coords)
                 or just a single atom.
                 which gives the list of atoms and the list of their coordinates
    basis_set    String describing the gaussian basis set to use
    ordering     The ordering of the basis functions to use
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

  def evaluate_at(self, x, y, z, mask=None):
    raise NotImplementedError("To be implemented")
    pass

