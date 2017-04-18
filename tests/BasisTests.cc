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
#include <gint/Structure.hh>
#include <gint/gaussian/Basis.hh>
#include <gint/gaussian/BasisSet.hh>
#include <gint/gaussian/Shell.hh>

namespace gint {
namespace tests {
using namespace gaussian;

TEST_CASE("basis tests", "[Basis]") {
  const BasisSet sto3g = lookup_basisset("sto-3g");
  const BasisSet g631 = lookup_basisset("6-31G*");

  const Structure h2{
        {"H", {{0, 0, 0}}},  //
        {"H", {{0, 0, 2}}},
  };

  const Structure h2o{
        {"O", {{0, 0, 0}}},
        {"H", {{0, 0, 1.795239827225189}}},
        {"H", {{1.693194615993441, 0, -0.599043184453037}}},
  };

  SECTION("Test Hydrogen molecule in sto3g") {

    // Construct a basis for h2 using the sto3g basis set:
    const Basis bas(h2, sto3g);
    CHECK(bas.size() == 2);

    // Check coords:
    for (size_t i = 0; i < h2.size(); ++i) {
      for (size_t j = 0; j < 3; ++j) {
        CHECK(bas[i].origin[j] == h2[i].coords[j]);
      }
    }

    // Check coefficients and exponents
    for (size_t i = 0; i < 2; ++i) {
      CHECK(bas[i].coefficients.size() == 3);
      CHECK(bas[i].exponents.size() == 3);

      CHECK(bas[i].coefficients[0] == 0.15432897);
      CHECK(bas[i].coefficients[1] == 0.53532814);
      CHECK(bas[i].coefficients[2] == 0.44463454);
      CHECK(bas[i].exponents[0] == 3.42525091);
      CHECK(bas[i].exponents[1] == 0.62391373);
      CHECK(bas[i].exponents[2] == 0.16885540);
    }
  }  // H2 sto-3g

  SECTION("Test Hydrogen molecule in 6-31G*") {
    const Basis bas(h2, g631);
    CHECK(bas.size() == 4);

    // Check coefficients and exponents
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        CHECK(bas[2 * i].origin[j] == h2[i].coords[j]);
      }

      CHECK(bas[2 * i].coefficients.size() == 3);
      CHECK(bas[2 * i].coefficients[0] == 0.03349460);
      CHECK(bas[2 * i].coefficients[1] == 0.23472695);
      CHECK(bas[2 * i].coefficients[2] == 0.81375733);

      CHECK(bas[2 * i].exponents.size() == 3);
      CHECK(bas[2 * i].exponents[0] == 18.7311370);
      CHECK(bas[2 * i].exponents[1] == 2.8253937);
      CHECK(bas[2 * i].exponents[2] == 0.6401217);

      for (size_t j = 0; j < 3; ++j) {
        CHECK(bas[2 * i + 1].origin[j] == h2[i].coords[j]);
      }

      CHECK(bas[2 * i + 1].coefficients.size() == 1);
      CHECK(bas[2 * i + 1].coefficients[0] == 1.0000000);

      CHECK(bas[2 * i + 1].exponents.size() == 1);
      CHECK(bas[2 * i + 1].exponents[0] == 0.1612778);
    }
  }  // H2 6-31G*

  SECTION("Test water 6-31G*") {
    const Basis bas(h2o, g631);
    CHECK(bas.size() == 10);

    // Oxygen
    for (size_t i = 0; i < 6; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        CHECK(bas[i].origin[j] == h2o[0].coords[j]);
      }
    }

    CHECK(bas[0].exponents.size() == 6);
    CHECK(bas[0].coefficients.size() == 6);
    CHECK(bas[0].l == 0);

    CHECK(bas[1].exponents.size() == 3);
    CHECK(bas[1].coefficients.size() == 3);
    CHECK(bas[1].l == 0);

    CHECK(bas[2].exponents.size() == 3);
    CHECK(bas[2].coefficients.size() == 3);
    CHECK(bas[2].l == 1);

    CHECK(bas[3].exponents.size() == 1);
    CHECK(bas[3].coefficients.size() == 1);
    CHECK(bas[3].l == 0);

    CHECK(bas[4].exponents.size() == 1);
    CHECK(bas[4].coefficients.size() == 1);
    CHECK(bas[4].l == 1);

    CHECK(bas[5].exponents.size() == 1);
    CHECK(bas[5].coefficients.size() == 1);
    CHECK(bas[5].l == 2);

    // Hydrogen 1
    CHECK(bas[6].exponents.size() == 3);
    CHECK(bas[6].coefficients.size() == 3);
    CHECK(bas[6].l == 0);
    CHECK(bas[7].exponents.size() == 1);
    CHECK(bas[7].coefficients.size() == 1);
    CHECK(bas[7].l == 0);

    for (size_t j = 0; j < 3; ++j) {
      CHECK(bas[6].origin[j] == h2o[1].coords[j]);
      CHECK(bas[7].origin[j] == h2o[1].coords[j]);
    }

    // Hydrogen 2
    CHECK(bas[8].exponents.size() == 3);
    CHECK(bas[8].coefficients.size() == 3);
    CHECK(bas[8].l == 0);
    CHECK(bas[9].exponents.size() == 1);
    CHECK(bas[9].coefficients.size() == 1);
    CHECK(bas[9].l == 0);
    for (size_t j = 0; j < 3; ++j) {
      CHECK(bas[8].origin[j] == h2o[2].coords[j]);
      CHECK(bas[9].origin[j] == h2o[2].coords[j]);
    }
  }  // H2O 6-31G*
}  // Basis

}  // namespace tests
}  // namespace gint
