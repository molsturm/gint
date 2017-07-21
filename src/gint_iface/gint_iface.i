//
// Copyright (C) 2017 by the gint authors
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

%module gint_iface

%{
#include <krims/ExceptionSystem.hh>

#include <gint/gaussian/Shell.hh>
#include <gint/gaussian/Basis.hh>
#include "construct_gaussian_basis.hh"

// Run the %init block (to setup numpy)
#define SWIG_FILE_WITH_INIT
%}

// TODO Extremely rudimentary exception handling
%exception {
  try {
    $action
  } catch (const krims::ExcNotImplemented& e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    return NULL;
  } catch (const krims::ExceptionBase& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }
}

%include "numpy.i"
%include "std_string.i"
%include "flat_structure.i"

// Setup import of numpy array:
%init %{
import_array();
%}

%include "../gint/gaussian/Shell.hh"
%include "../gint/gaussian/Basis.hh"
%include "construct_gaussian_basis.hh"

// vi: syntax=c
