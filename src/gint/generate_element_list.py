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

from datetime import date
import os
import sys


# Amend path to make sure that the element.py file is actually found.
sys.path = ["../gint_iface/"] + sys.path
from gint.element import elements


OUTPUT_FILE = "elements.cc"
HEADER = """//
// Copyright (C) """ + str(date.today().year) + """ by the gint authors
//
// This file is part of gint.
//
// gint is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gint is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gint. If not, see <http://www.gnu.org/licenses/>.
//

//
// Note: This file is generated by the script """ + os.path.basename(__file__) + """
//       Do not modify this file directly, but edit the generator
//       script instead.
//

#include "elements.hh"

namespace gint {

const std::vector<Element>& elements() {
  static std::vector<Element> elements{{
"""

FOOTER = """
  }};

  return elements;
}

}  // namespace gint
"""

# Build format string for the individual element lines
maxlen = max(len(e.name) for e in elements) + 2
fmt = '        {{{atom_number:3d}, {symbol:5s}, {name:' + str(maxlen) + 's}}},\n'

# Build the element lines
elementlines = []
for element in elements:
    params = {
        "name":        '"' + element.name + '"',
        "symbol":      '"' + element.symbol + '"',
        "atom_number": element.atom_number,
    }
    elementlines.append(fmt.format(**params))

# Strip the last newline
elementlines[-1] = elementlines[-1][:-1]

# Dump the content
with open(OUTPUT_FILE, "w") as f:
    f.write(HEADER)
    f.writelines(elementlines)
    f.write(FOOTER)
