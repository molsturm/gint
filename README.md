# gint
[![Build Status](https://travis-ci.org/molsturm/gint.svg?branch=master)](https://travis-ci.org/molsturm/gint)
[![Coverage Status](https://coveralls.io/repos/github/molsturm/gint/badge.svg?branch=master)](https://coveralls.io/github/molsturm/gint)
[![Licence](https://img.shields.io/github/license/molsturm/gint.svg)](LICENCE)

A basis-function independent interface library providing the types of
integrals needed for electronic structure theory calculations.

## Design goals
- The same interface independent of the type of basis function or the precise
  integral backend.

More documentation will follow ...

## Building and interfacing with external libraries
Currently `gint` can be build with a couple of external libraries for performing
the integral computation

### libint
Edward Valeev's [``libint``](https://github.com/evaleev/libint) library
can be enabled as the Gaussian integral backend by passing
`-DGINT_ENABLE_LIBINT=ON` to `cmake` at configure time,
which will automatically checkout and build `libint` along with `gint`.

Further external dependencies required by `libint`:
- Eigen3
- BLAS
- Autoconf
- GNU Multiprecision library

On Debian/Ubuntu the command
```
apt-get install libeigen3-dev libblas-dev autoconf libgmp-dev
```
should do.

### libcint
The second supported Gaussian integral library is [``libcint``](https://github.com/sunqm/libcint)
by Qiming Sun. It can be enabled by `-DGINT_ENABLE_LIBCINT=ON`,
which will again compile `libcint` along `gint`.

Further external dependencies required by `libcint`:
- BLAS

### sturmint
Sturmian and Coulomb-Sturmian library. Currently not publicly available (but planned to be released soon).
Can be enabled by `-DGINT_ENABLE_STURMINT=ON`.
