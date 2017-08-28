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

#include "available_basis_types.hh"
#include <gint/IntegralLookup.hh>
#include <gint/config.hh>

namespace gint {
namespace interface {

std::vector<std::string> available_basis_types() {
  using real_type = gint::real_valued::stored_matrix_type;
  using cpx_type  = gint::complex_valued::stored_matrix_type;

  // Combine complex and real basis function types.
  std::vector<std::string> ret = gint::IntegralLookup<real_type>::available_basis_types();
  const std::vector<std::string> cpx =
        gint::IntegralLookup<cpx_type>::available_basis_types();
  std::copy(cpx.begin(), cpx.end(), std::back_inserter(ret));

  // Sort and unique
  std::sort(ret.begin(), ret.end());
  std::unique(ret.begin(), ret.end());
  return ret;
}

}  // namespace interface
}  // namespace gint
