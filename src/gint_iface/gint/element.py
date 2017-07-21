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

import collections
import collections.abc

Element = collections.namedtuple("Element", [ "atom_number", "symbol", "name" ])

# TODO Autogenerate from shared data with the C++ side
elements = (
  Element(  1, "H",  "hydrogen"  ),
  Element(  2, "He", "helium"    ),
  Element(  3, "Li", "lithium"   ),
  Element(  4, "Be", "beryllium" ),
  Element(  5, "B",  "boron"     ),
  Element(  6, "C",  "carbon"    ),
  Element(  7, "N",  "nitrogen"  ),
  Element(  8, "O",  "oxygen"    ),
  Element(  9, "F",  "fluorine"  ),
  Element( 10, "Ne", "neon"      ),
  Element( 11, "Na", "sodium"    ),
  Element( 12, "Mg", "magnesium" ),
  Element( 13, "Al", "aluminium" ),
  Element( 14, "Si", "silicon"   ),
  Element( 15, "P",  "phosphorus"),
  Element( 16, "S",  "sulphur"   ),
  Element( 17, "Cl", "chlorine"  ),
  Element( 18, "Ar", "argon"     ),
)


def by_symbol(symbol):
  """
  Lookup an element by element symbol.
  Returns the Element datastructure of the object found.

  Raises a KeyError if the element cannot be found.
  """
  symbol = symbol.lower()
  el = [ el for el in elements if el.symbol.lower() == symbol ]
  if len(el) != 1:
    raise KeyError("No element with the symbol " + symbol + " found.")
  return el[0]


def by_atomic_number(num):
  """
  Lookup an element by atomic number.
  Returns the Element datastructure of the object found.

  Raises a KeyError if the element cannot be found.
  """
  el = [ el for el in elements if el.atom_number == num ]
  if len(el) != 1:
    raise KeyError("No element with the atomic number " + str(num) + " found.")
  return el[0]


def by_name(name):
  """
  Lookup an element by name. Returns the Element
  datastructure of the object found.

  Raises a KeyError if the element cannot be found.
  """
  name = name.lower()
  el = [ el for el in elements if el.name.lower() == name ]
  if len(el) != 1:
    raise KeyError("No element with the name " + name + " found.")
  return el[0]


def find(val):
  """
  Lookup an element by some criterion. The criterion
  is interpreted as follows:

  - If it is an integer, by_atomic_number is called
  - If it is a string with 3 or less characters,
    it is interpreted as an element symbol
  - Else it is interpreted as an element name.
  """
  if isinstance(val, int):
    return by_atomic_number(val)
  elif isinstance(val, str):
    if len(val) > 3:
      return by_name(val)
    else:
      return by_symbol(val)
  else:
    raise TypeError("val should be an integer or string")

by_anything = find


def to_atom_numbers(elements):
  """
  Flexible function for converting to element numbers.
  It essentially calls ``by_anything(item).to_atom_number``
  on each of the items of the iterable.
  """
  if isinstance(elements, collections.abc.Iterable):
    return [ by_anything(item).atom_number for item in elements ]
  else:
    return [ by_anything(elements).atom_number ]


