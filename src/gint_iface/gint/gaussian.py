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

Shell = collections.namedtuple("Shell", [ "l", "pure", "contraction_coefficients",
                                          "contraction_exponents", "centre" ])


class Basis:
  def __init__(self, atoms, coords, basis_set):
    if len(atoms) != len(coords):
      raise ValueError("Number of atoms and number of coordinates does not agree.")

    atom_numbers = np.array(element.to_atom_numbers(atoms))
    coords = np.array(coords)
    if coords.shape[1] != 3:
      raise ValueError("The coords list needs to have exactly 3 items per coordinate.")

    basis_raw = _iface.construct_gaussian_basis(atom_numbers, coords, basis_set)
    shells_raw = ( basis_raw.shell(i) for i in range(basis_raw.n_shells()) )

    # Setup the shells.
    # Note the explicit copy! This is needed, since the basis_raw.shell
    # function merely returns a view into the data, which is held by
    # the c++ object behing basis_raw.
    self.shells = [ Shell(sh.l, sh.pure, sh.coefficients().copy(),
                          sh.exponents().copy(), sh.origin().copy())
                    for sh in shells_raw ]

  def evaluate_at(self, x, y, z, mask=None):
    pass

