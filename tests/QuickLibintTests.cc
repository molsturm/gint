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

#ifdef GINT_HAVE_LIBINT

#include "GaussianTestData.hh"
#include "gint/config.hh"
#include "integral_quick_tests.hh"
#include <catch.hpp>
#include <gint/IntegralLookup.hh>
#include <linalgwrap/SmallVector.hh>
#include <rapidcheck.h>

namespace gint {

namespace tests {
using namespace linalgwrap;
using namespace krims;

TEST_CASE("Quick atomic libint test", "[quicktest libint]") {
  const OrbitalType otype = OrbitalType::REAL_MOLECULAR;
  typedef IntegralLookup<otype> int_lookup_type;

  // The reference data for atomic coulomb sturmians
  // with parameters k = 1, Z = 4, n_max =  3, l_max = 2
  typedef GaussianTestData<real_valued::stored_matrix_type> refdata_type;

  // TODO This is bad, but to get going ....
  auto highertol = krims::NumCompConstants::change_temporary(
        1e6 * krims::NumCompConstants::default_tolerance_factor);

  // Setup parameters for the integral library
  const krims::GenMap params{
        {"basis_type", "gaussian/libint"},
        {"basis_set", refdata_type::basis},
        {"structure", refdata_type::molecule},
  };

  IntegralDummyTests<int_lookup_type, refdata_type>::run_all("libint2: ",
                                                             int_lookup_type(params));
}

}  // namespace tests
}  // namespace gint

#endif  // GINT_HAVE_LIBINT
