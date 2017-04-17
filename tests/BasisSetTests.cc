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
#include <gint/gaussian/BasisSet.hh>
#include <gint/gaussian/Shell.hh>
#include <krims/FileUtils.hh>

namespace gint {
namespace tests {
using namespace gaussian;

TEST_CASE("BasisSet tests", "[basis_set]") {
  SECTION("Test lookup_basisset for sto-3g") {
    const BasisSet bs = lookup_basisset("sto-3g");

    CHECK(bs.name == "sto-3g");
    CHECK(bs.filename.substr(bs.filename.size() - 16) == "basis/sto-3g.g94");
    CHECK(bs.atomic_number_to_shells.size() == 53);

    SECTION("Test hydrogen") {
      auto& h = bs.shells_for_atom(1);
      CHECK(h.size() == 1);
      CHECK(h[0].coefficients[0] == 0.15432897);
      CHECK(h[0].coefficients[1] == 0.53532814);
      CHECK(h[0].coefficients[2] == 0.44463454);
      CHECK(h[0].exponents[0] == 3.42525091);
      CHECK(h[0].exponents[1] == 0.62391373);
      CHECK(h[0].exponents[2] == 0.16885540);
    }

    // Magnesium
    CHECK(bs.shells_for_atom(12).size() == 5);

    // Silver
    CHECK(bs.shells_for_atom(47).size() == 11);

    // Osmium and Xenon (should not exist)
    CHECK_THROWS_AS(bs.shells_for_atom(54), ExcNoBasisForAtom);
    CHECK_THROWS_AS(bs.shells_for_atom(76), ExcNoBasisForAtom);
  }  // sto-3g

  SECTION("Test lookup_basisset for 6-31g*") {
    const BasisSet bs = lookup_basisset("6-31g*");

    CHECK(bs.name == "6-31g*");
    CHECK(bs.filename.substr(bs.filename.size() - 16) == "basis/6-31g*.g94");
    CHECK(bs.atomic_number_to_shells.size() == 30);

    // Hydrogen
    CHECK(bs.shells_for_atom(1).size() == 2);

    // Magnesium
    CHECK(bs.shells_for_atom(12).size() == 8);

    // Cobalt
    CHECK(bs.shells_for_atom(27).size() == 12);
    CHECK(bs.shells_for_atom(27).back().l == 3);  // Check for F Shell

    // Osmium and Xenon (should not exist)
    CHECK_THROWS_AS(bs.shells_for_atom(54), ExcNoBasisForAtom);
    CHECK_THROWS_AS(bs.shells_for_atom(76), ExcNoBasisForAtom);
  }  // sto-3g

  SECTION("Check normalisation") {
    const BasisSet bs1 = lookup_basisset("sto-3g");
    const BasisSet bs2 = lookup_basisset("sto-3G");
    const BasisSet bs3 = lookup_basisset("sTo-3G");
    CHECK(bs1.name == bs2.name);
    CHECK(bs1.name == bs3.name);
  }

  SECTION("Check for failing examples") {
    CHECK_THROWS_AS(lookup_basisset("iaen"), krims::ExcDatafileNotFound);
  }

  // TODO Test augmentation basis sets like aug-cc-pvdz
  //      (which are read from 2 files !)

}  // BasisSet

}  // namespace tests
}  // namespace gint
