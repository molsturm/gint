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
import numpy as np

class Structure:
  def __init__(self, atoms, coords=None):
    """Setup a structure object.

    atoms:   List of all atoms of the structure. Either symbol or atomic numbers
             or names are accepted.
    coords:  List of the coordinates of the involved atoms.

    If there is only a single atom, then coords may be absent and atoms may be just
    a single atom number or atom name/symbol.
    """
    if not isinstance(atoms, (list, tuple)):
      if coords is None:
        atoms = [ atoms ]
      else:
        raise TypeError("atoms needs to be a list or tuple")

    if coords is None:
      if len(atoms) == 1:
        coords = [[0,0,0]]
      else:
        raise ValueError("coords needs to be specified if more than one atom is in "
                         "the structure.")

    if len(atoms) != len(coords):
      raise ValueError("Number of atoms and number of coordinates does not agree.")

    self.atom_numbers = np.array(element.to_atom_numbers(atoms))
    self.coords = np.array(coords)
    if self.coords.shape[1] != 3:
      raise ValueError("The coords list needs to have exactly 3 items per coordinate.")


  @property
  def atoms(self):
    """Return the list of elements which make up the structure"""
    return [ element.by_atomic_number(i) for i in self.atom_numbers ]


  @property
  def n_atoms(self):
    return len(self.atom_numbers)


  @property
  def total_charge(self):
    """Compute the total charge of the atoms in the structure"""
    return sum(self.atom_numbers)


