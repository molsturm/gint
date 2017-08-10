#!/usr/bin/env python3
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab
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

from datetime import date
import os
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
    "φ":    "phi",
    "θ":    "theta",
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

  For this normalisation factor ζ is arbitrary and included
  in the expression.
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


def normalisation_pure_gaussians(l):
  """
  Generates the factor which normalises the radial part of a
  Gaussian with spherial harmonics, i.e. such that
  r**2 * norm**2 * func**2 will integrate to 1 over r = (0, oo)
  where func = exp(-r*r)

  Note that here we do *not* include the ζ, i.e. we implicitly
  assume ζ == 1 here.
  For a ζ != 1, the returned expression needs to be scaled by
  ζ**((2*l+3)/4) in order to be correct.
  """
  r = sy.Symbol("r", positive=True)

  f = r**l * sy.exp(-r*r)
  norm = sy.integrate( f*f*r*r, (r,0, sy.oo) )
  return 1 / sy.sqrt(norm)


def spherical_harmonics_times_r(l):
  """
  Generates the terms for r**l * Y_lm(l, m, θ, φ) for m in range(-l, l+1)
  and returns them as a list. Y_lm in this context are *real* spherical
  harmonics.

  The returned terms T are normalised such that
  ( T * ζ**((2*l+3)/4) * exp(- ζ*r*r) )**2
  integrates to unity over space.
  """
  r = sy.Symbol("r", positive=True)
  θ = sy.Symbol("θ", real=True)
  φ = sy.Symbol("φ", positive=True)
  x, y, z = sy.symbols("x, y, z", real=True)

  # Cartesian substitutions
  sub_for_cart = {
    r * sy.sin(θ) * sy.cos(φ) : x,
    r * sy.sin(θ) * sy.sin(φ) : y,
    r * sy.cos(θ)             : z,
  }

  # The normalisation for this shell
  norm = normalisation_pure_gaussians(l)

  ret = []
  for m in range(-l, l+1):
    rYlm = norm * sy.expand(r**l * sy.Znm(l, m, θ, φ), func=True)
    rYlm, imag = rYlm.as_real_imag()
    assert imag == 0
    rYlm = sy.simplify(rYlm)
    rYlm = rYlm.subs(sub_for_cart)
    ret.append(rYlm)
  return ret

#
# -----------------------------------------------------------------
#

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

def code_angular_pack_pure_gaussian(lmax):
  string = "  angular=None\n"
  for l in range(lmax+1):
    el = "el" if l > 0 else ""
    string += "  " + el  + "if l == " + str(l) + ":\n"
    string += "    angular=[\n"
    for mexpr in spherical_harmonics_times_r(l):
      string += "        " + as_numpy_code(mexpr) + ",\n"
    string += "    ]\n"
  string += "  else:\n"
  string += "    raise NotImplementedError(\"Only pure gaussians up to l=="
  string += str(lmax) + " are implemented.\")\n"
  return string

#
# -----------------------------------------------------------------
#

def dump_shell_expression_file(file, l_max):
  """
  Dump the python file holding the expression needed for evaluating
  shells of gaussian functions up to a given l_max.
  """

  with open(file, "w") as f:
    f.write("## "+69*"-"+"\n")
    f.write("##\n")
    f.write("## Copyright (C) " +str(date.today().year) + " the gint authors\n")
    f.write("##\n")
    f.write("## This file is licenced under the GNU General Public License\n")
    f.write("## version 3. See the other files the gint project for more details.\n")
    f.write("##\n")
    f.write("## This file is automatically generated. Do not edit directly.\n")
    f.write("## Instead edit the generator script \"" + os.path.basename(__file__) + "\"\n")
    f.write("##\n")
    f.write("## "+69*"-"+"\n\n")

    f.write("import numpy as np\n\n")

    f.write("def normalisation_cartesian_gaussians(l, zeta):\n")
    f.write(code_normalisation_cartesian_gaussian(l_max))
    f.write("  return norm\n")
    f.write("\n\n")

    f.write("def angular_pack_pure_gaussian(l, x, y, z, r, theta, phi):\n")
    f.write(code_angular_pack_pure_gaussian(l_max))
    f.write("  return angular\n")
    f.write("\n\n")


if __name__ == "__main__":
  dump_shell_expression_file("_gaussian_shells_expressions.py", l_max=7)

