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

## Dependencies
For *building* `gint` the following things are required:
- ``cmake`` >= 3.0.0
- A compiler supporting ``C++11``: ``clang`` starting from `clang-3.5` and `gcc` starting
  from `gcc-4.8` should work.
- [``swig``](http://swig.org/) >= 2.0.11
- [``python``](https://www.python.org/) >= 3.4, including the development headers
- The [usual build process](#building-gint) mentioned below
  will automatically build the [``lazyten``](https://lazyten.org) linear algebra library
  as well. This requires further
    - A BLAS implementation, e.g. [OpenBLAS](https://github.com/xianyi/OpenBLAS/)
    - A LAPACK compatible library, e.g. [LAPACK](http://netlib.org/lapack)
    - [armadillo](http://arma.sourceforge.net/)

  See [github.com/lazyten/lazyten](https://github.com/lazyten/lazyten/blob/master/README.md)
  for more details about ``lazyten``'s dependencies.
- If `gint` should be linked with third-party libraries to compute any integrals,
  the requirements of these will be needed as well.
  See the section [Interfacing with external libraries](#interfacing-with-external-libraries).

In order to actually *use* the `gint` python module once it has been built,
the following `python` packages are required:
- [numpy](https://pypi.python.org/pypi/numpy)


On a recent **Debian/Ubuntu** you can install all these dependencies by running
```
apt-get install cmake swig python3-dev libopenblas-dev liblapack-dev libarmadillo-dev \
                python3-numpy
```
as root.

## Building ``gint``
This basic build builds `gint` along with one or more
third-party libraries to compute the integrals as well as tests.
The basic sequence of steps is
```sh
# Configure inside build directory
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=~/opt -DAUTOCHECKOUT_MISSING_REPOS=ON ${OPTIONS} ..
cmake --build .

# Test and install the build
ctest
cmake --build . --target install
```
where `${OPTIONS}` is one of the `cmake` options mentioned in the next section.

## Interfacing with external libraries
Currently `gint` can be build with a couple of external libraries for performing
the integral computation.

### Gaussians: libint
Edward Valeev's [``libint``](https://github.com/evaleev/libint) library
is enabled by passing the option `-DGINT_ENABLE_LIBINT=ON` to `cmake`.
This will automaticall download, unpack, compile and link to `libint`.

Next to the dependencies mentioned before, building `libint` requires:
- Eigen3
- BLAS
- Autoconf
- GNU Multiprecision library

On Debian/Ubuntu the command
```
apt-get install libeigen3-dev libblas-dev autoconf libgmp-dev
```
should install all of these.

### Gaussians: libcint
The second supported Gaussian integral library is [``libcint``](https://github.com/sunqm/libcint)
by Qiming Sun. It can be enabled by the option `-DGINT_ENABLE_LIBCINT=ON`,
which will again compile `libcint` along `gint`.

Further external dependencies required by `libcint`:
- BLAS

### Sturmians: sturmint
Sturmint is currently not released (but this is planned for the near future).
It can be enabled using the option `-DGINT_ENABLE_STURMINT=ON`.
