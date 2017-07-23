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

import numpy as np
from scipy.misc import factorial2

def cartesian_xyz_expontents(l, ordering):
  """
  Produce the list of cartesian exponents.
  This determines the order inside the cartesian exponents, which is used.

  The precise ordering can be adapted using the ordering parameter.
  The following orderings are currently implemented:

  - standard: The ordering of the Common Component Architecture (CCA) standard
              as described in the paper DOI 10.1002/jcc
  """
  # Special cases => deal with those first
  if l == 0:
    return [ (0, 0, 0) ]
  elif l == 1:
    return [ (1, 0, 0), (0, 1, 0), (0, 0, 1) ]

  if ordering == "standard":
    # Standard ordering as described in DOI 10.1002/jcc
    return [ (a, l-a-c, c) for a in range(l, -1, -1) for c in range(0, l-a+1) ]
  else:
    raise NotImplementedError("No ordering except standard is implemented.")


def cartesian_shell(x, y, z, l, centre, coefficient, zeta, ordering):
  """
  Return the function values of a shell of contracted cartesian
  gaussians with orbital angular momentum l and exponent zeta
  centered at centre.

  @param x       The x values on the mesh
  @param y       The y values on the mesh
  @param z       The z values on the mesh
  @param l       The angular momentum of the shell
  @param centre  The centre of the Gaussians
  @param coefficient  np array of contraction coefficients
  @param zeta         np array of contraction exponents
  @param ordering  The ordering to use inside the cartesian harmonics.
  """
  if len(centre) != 3:
    raise ValueError("centre needs to have exactly 3 components")
  if len(coefficient) != len(zeta):
    raise ValueError("Length of contraction coefficients and length of contraction "
                     "exponents(zeta) needs to agree")


  # According to DOI 10.1002/jcc norm of a primitive is the factor
  # normalising the cartesian function
  #     z^l exp(-ζ * r*r) = (r cos θ)^l exp(-ζ * r*r)
  # in three-dimensional space. For various ranges of l presented
  # here the expressions have been computed using
  # generate_gaussian_expressions.py
  # -------------------------------------------------------
  norm=None
  if l == 0:
    norm = 2**(3/4)*zeta**(3/4)/np.pi**(3/4)
  elif l == 1:
    norm = 2*2**(3/4)*zeta**(5/4)/np.pi**(3/4)
  elif l == 2:
    norm = 4*2**(3/4)*np.sqrt(3)*zeta**(7/4)/(3*np.pi**(3/4))
  elif l == 3:
    norm = 8*np.sqrt(15)*2**(3/4)*zeta**(9/4)/(15*np.pi**(3/4))
  elif l == 4:
    norm = 16*np.sqrt(105)*2**(3/4)*zeta**(11/4)/(105*np.pi**(3/4))
  elif l == 5:
    norm = 32*np.sqrt(105)*2**(3/4)*zeta**(13/4)/(315*np.pi**(3/4))
  elif l == 6:
    norm = 64*np.sqrt(1155)*2**(3/4)*zeta**(15/4)/(3465*np.pi**(3/4))
  elif l == 7:
    norm = 128*np.sqrt(15015)*2**(3/4)*zeta**(17/4)/(45045*np.pi**(3/4))
  elif l == 8:
    norm = 256*np.sqrt(1001)*2**(3/4)*zeta**(19/4)/(45045*np.pi**(3/4))
  elif l == 9:
    norm = 512*np.sqrt(17017)*2**(3/4)*zeta**(21/4)/(765765*np.pi**(3/4))
  elif l == 10:
    norm = 1024*2**(3/4)*np.sqrt(323323)*zeta**(23/4)/(14549535*np.pi**(3/4))
  else:
    raise NotImplementedError("Only cartesian gaussians up to l==10 are implemented.")
  # -------------------------------------------------------

  if np.all(centre != np.zeros(3)):
    xs = x - centre[0]
    ys = y - centre[1]
    zs = z - centre[2]
  else:
    (xs, ys, zs) = (x, y, z)

  # Compute the radial part by broadcasting the norms, the zetas and the xs, ys, zs
  # such that the computation can occurr in one step. Contract the resulting tensor
  # on the first axis thereafter
  prefac = norm * coefficient

  prefac = prefac[:, np.newaxis, np.newaxis, np.newaxis]
  zeta   = zeta[:, np.newaxis, np.newaxis, np.newaxis]
  xsb    = xs[np.newaxis, :, :, :]
  ysb    = ys[np.newaxis, :, :, :]
  zsb    = zs[np.newaxis, :, :, :]
  radial = prefac * np.exp(- zeta * (xsb**2 + ysb**2 + zsb**2))
  radial = np.sum(radial, axis=0)

  # Now compute the angular part.
  # For this keep in mind that we need to adjust the normalisation for each
  # function individually
  xyz_exponents = cartesian_xyz_expontents(l, ordering=ordering)
  result = np.empty( (len(xyz_exponents),) + x.shape, dtype=x.dtype )
  for i, (px, py, pz) in enumerate(xyz_exponents):
    ang_norm = factorial2(2*px + 2*py + 2*pz - 1)
    ang_norm /= (factorial2(2*px - 1) * factorial2(2*py - 1) * factorial2(2*pz - 1))
    result[i, :, :, :] = ang_norm * xs**px * ys**py * zs**pz * radial
  return result


def pure_shell(x, y, z, l, centre, coefficient, zeta):
  raise NotImplementedError("Pure gaussians are not yet implemented")
