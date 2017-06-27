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

#include <catch.hpp>
#include <gint/IntegralLookup.hh>
#include <gint/config.hh>

namespace gint {
namespace tests {

TEST_CASE("IntegralLookup tests", "[integral_lookup]") {
  SECTION("Test available_basis_types()") {
    std::vector<std::string> avail =
          IntegralLookup<real_valued::stored_matrix_type>::available_basis_types();

    auto find_value = [&avail](const std::string& value) {
      return std::find(avail.begin(), avail.end(), value) != avail.end();
    };
    CHECK_FALSE(find_value("atomic/cs_naive"));

#ifdef GINT_HAVE_STATIC_INTEGRALS
    CHECK(find_value("sturmian/atomic/cs_static14"));
#endif  // GINT_HAVE_STATIC_INTEGRALS

#ifdef GINT_HAVE_STURMINT
    CHECK(find_value("sturmian/atomic/cs_reference"));
    CHECK(find_value("sturmian/atomic/cs_naive"));
#endif  // GINT_HAVE_STURMINT

#ifdef GINT_HAVE_LIBINT
    CHECK(find_value("gaussian/libint"));
#endif  // GINT_HAVE_LIBINT

  }  // Test available

  // TODO Trigger and test registering basis types
}  // Lookup tests

}  // namespace tests
}  // namespace gint
