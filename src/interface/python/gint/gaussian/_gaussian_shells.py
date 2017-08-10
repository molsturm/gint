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

import numpy as np
from scipy.misc import factorial2
from ..util import cartesian_to_spherical
from . import _gaussian_shells_expressions as expr


def norm_contraction(l, coefficient, zeta):
  """
  In Gaussian basis sets the contraction coefficients are usually
  given relative to normalised primitive functions.

  Since the primitives are not orthogonal, however, the overlap
  between them needs to be taken into account when forming
  the contraction by linear combination.

  Some basis sets do this by incorating this into the values of the
  contraction coefficients, but some don't.

  This function computes the norm factor required. Assume that
  φ is a contracted gaussian, i.e.
       φ = \sum_i c_i f_i / \sqrt{<f_i | f_i>}
  where c_i are the contraction coefficients and f_i are the primitves.
  Then
       <φ|φ> = \sum_{ij} c_i c_j <f_i | f_j> / \sqrt{<f_i|f_i> <f_j|f_j>}
  Since f_i = r^l * exp(-ζ_i * r^2) we now find
       ovlap_{ij} := <f_i | f_j> / \sqrt{<f_i|f_i> <f_j|f_j>}
                   = (4 ζ_i ζ_j / (ζ_i + ζ_j)^2) ^ ((2*l + 3)/4)
  In other words the required norm factor is given as
       1 / sqrt(<φ|φ>)
  where <φ|φ> is the contraction-coefficient weighted sum of
  the overlap factors ovlap_{ij} introduced above.
  """
  zeta_i = zeta[:, np.newaxis]
  zeta_j = zeta[np.newaxis, :]
  coeff_i = coefficient[:, np.newaxis]
  coeff_j = coefficient[np.newaxis, :]

  lpow = (2 * l + 3) / 4
  ovlap_ij = (4 * zeta_i * zeta_j / (zeta_i + zeta_j) ** 2) ** lpow
  return 1 / np.sqrt(np.sum(coeff_i * coeff_j * ovlap_ij))


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
  if centre.ndim != 1:
    raise ValueError("centre needs to be a 1-D numpy array.")

  if len(coefficient) != len(zeta):
    raise ValueError("Length of contraction coefficients and length of contraction "
                     "exponents(zeta) needs to agree")
  if coefficient.ndim != 1 or zeta.ndim != 1:
    raise ValueError("coeffients and zeta need to be 1-D numpy arrays.")

  # According to DOI 10.1002/jcc norm of a primitive is the factor
  # normalising the cartesian function
  #     z^l exp(-ζ * r*r) = (r cos θ)^l exp(-ζ * r*r)         (1)
  # in three-dimensional space. For various ranges of l presented
  # here the expressions have been computed using the script
  # _generate_shell_expressions.py and put into the output
  # _shell_expressions.py
  norm = expr.normalisation_cartesian_gaussians(l, zeta)

  if np.all(centre != np.zeros(3)):
    xs = x - centre[0]
    ys = y - centre[1]
    zs = z - centre[2]
  else:
    (xs, ys, zs) = (x, y, z)

  # In order to normalise the contracted gaussian we need
  #  - The contraction factor, which is the norm change due to the contraction
  #    relative to *normalised* primitives
  #  - The norm factor which actually normalise the primitives by themselves
  #    (i.e. the thing computed above together with an adjusting factor
  #    ang_norm --- see below, that adjusts the norm of a particular Cartesian
  #    harmonic relative to the norm of the standard function given in (1) above.)
  prefac = norm * coefficient * norm_contraction(l, coefficient, zeta)

  # Compute the radial part by broadcasting the norms, the zetas and the xs, ys, zs
  # such that the computation can occur in one step. We sum the resulting tensor
  # on the first axis thereafter to build the contracted gaussian linear combination.
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
    ang_norm = np.sqrt(ang_norm)

    result[i, :, :, :] = ang_norm * xs**px * ys**py * zs**pz * radial
  return result


def pure_shell(x, y, z, l, centre, coefficient, zeta):
  """
  Return the function values of a shell of pure cartesian
  gaussians with orbital angular momentum l and exponent zeta
  centered at centre.

  A pure gaussian is a gaussian radial part combined with a spherical
  harmonic as an angular part.

  @param x       The x values on the mesh
  @param y       The y values on the mesh
  @param z       The z values on the mesh
  @param l       The angular momentum of the shell
  @param centre  The centre of the Gaussians
  @param coefficient  np array of contraction coefficients
  @param zeta         np array of contraction exponents
  """
  if len(centre) != 3:
    raise ValueError("centre needs to have exactly 3 components")
  if centre.ndim != 1:
    raise ValueError("centre needs to be a 1-D numpy array.")

  if len(coefficient) != len(zeta):
    raise ValueError("Length of contraction coefficients and length of contraction "
                     "exponents(zeta) needs to agree")
  if coefficient.ndim != 1 or zeta.ndim != 1:
    raise ValueError("coeffients and zeta need to be 1-D numpy arrays.")

  if np.min(zeta) <= 0:
    raise ValueError("All exponents (zeta) need to be larger than zero")

  if np.all(centre != np.zeros(3)):
    xs = x - centre[0]
    ys = y - centre[1]
    zs = z - centre[2]
  else:
    (xs, ys, zs) = (x, y, z)

  # Compute the radial coordinates corresponding to x, y, z
  if l > 1:
    r, theta, phi = cartesian_to_spherical(xs, ys, zs)
  else:
    r = None
    theta = None
    phi = None

  # In order to normalise the contracted gaussian we need
  #  - The contraction factor, which is the norm change due to the contraction
  #    relative to *normalised* primitives
  #  - The norm factor which actually normalise the primitives by themselves
  #
  # The latter norm is further split up such that all factors which are
  # independent of zeta will be evaluated together with the angular part.
  # Only the zeta-dependent factor ζ**( (2*l + 3)/4 ) is evaluated here
  # and multiplied into the contraction coefficients.
  lpow = (2 * l + 3) / 4
  prefac = zeta**lpow * coefficient * norm_contraction(l, coefficient, zeta)
  assert prefac.shape == zeta.shape and prefac.shape == coefficient.shape

  # Again the procedure is analogous to the cartesian gaussians,
  # such that we broadcast for the radial part and sum over the axis
  # which contains the various contraction coefficients.
  # Note that the factor r**l the radial part usually holds is also
  # absorbed into the angular part (such that substitutions r*cos(θ) => z
  # become possible)
  prefac = prefac[:, np.newaxis, np.newaxis, np.newaxis]
  zeta   = zeta[:, np.newaxis, np.newaxis, np.newaxis]
  rs_sq  = (xs**2 + ys**2 + zs**2)[np.newaxis, :, :, :]
  radial = prefac * np.exp(-zeta * rs_sq)
  radial = np.sum(radial, axis=0)
  del prefac
  del rs_sq

  # The pack of angular integrals (in m=-l to m=l) we need to form the final
  # functions below. These contain the remaining factors of the normalisation
  # and the angular parts in the form of spherical harmonics, too
  # The expressions are generated by _gaussian_shells_expressions.py and are
  # stored in the _shell_expressions module
  angular = expr.angular_pack_pure_gaussian(l, x, y, z, r, theta, phi)
  result  = np.empty( (len(angular),) + x.shape, dtype=x.dtype )
  for i, ang in enumerate(angular):
    result[i, :, :, :] = ang * radial
  return result


