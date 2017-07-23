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

import sympy as sy


def as_numpy_code(inp):
  inp = sy.simplify(inp)
  s = str(inp)

  replacements = {
    "sin":  "np.sin",
    "pi":   "np.pi",
    "cos":  "np.cos",
    "exp":  "np.exp",
    "sqrt": "np.sqrt",
    "I":    "1j",
    "ζ":    "zeta",
  }
  for r in replacements:
    s = s.replace(r, replacements[r])
  return s


def normalisation_cartesian_gaussians(l):
  """
  Compute the normalisation coefficient for cartesian gaussians,
  following the convention introduced in DOI 10.1002/jcc, which
  means that one should take the largest angular momentum on
  one cartesian direction to compute the value
  """
  r = sy.Symbol("r", positive=True)
  θ = sy.Symbol("θ", positive=True)
  φ = sy.Symbol("φ", positive=True)
  ζ = sy.Symbol("ζ", positive=True)

  z = r*sy.cos(θ)
  f = z**l * sy.exp(-ζ*r*r)

  dV = r*r*sy.sin(θ)
  norm = sy.integrate( f*f*dV, (r,0, sy.oo), (θ, 0, sy.pi), (φ, -sy.pi, sy.pi))
  return 1 / sy.sqrt(norm)


def code_normalisation_cartesian_gaussian(lmax):
  string = "  norm=None\n"
  for l in range(lmax+1):
    el = "el" if l > 0 else ""
    string += "  " + el  + "if l == " + str(l) + ":\n"
    string += "    norm = " + as_numpy_code(normalisation_cartesian_gaussians(l)) + "\n"
  string += "  else:\n"
  string += "    raise NotImplementedError(\"Only cartesian gaussians up to l=="
  string += str(lmax) + " are implemented.\")\n"
  return string


if __name__ == "__main__":
  lmax = 10
  print("# Norms for cartesian gaussians z^l exp(-zeta*r*r):")
  print("#")
  print(code_normalisation_cartesian_gaussian(lmax))

