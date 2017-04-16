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
#include <gint/find_data_file.hh>
#include <gint/gaussian/BasisSet.hh>
#include <gint/gaussian/Shell.hh>
#include <gint/gaussian/read_basisset.hh>

namespace gint {
namespace tests {

TEST_CASE("read_basisset tests", "[read_basisset]") {

  SECTION("sto-3g test") {
    std::string file = find_data_file("basis/sto-3g.g94");
    std::ifstream f(file);
    const gaussian::BasisSet b =
          gaussian::read_basisset(f, gaussian::BasisSetFileFormat::Gaussian94);

    SECTION("Check hydrogen") {
      const std::vector<gaussian::Shell>& sh_h = b.shells_for_atom(1);
      CHECK(sh_h.size() == 1);

      const gaussian::Shell& h_1s = sh_h[0];
      CHECK(h_1s.l == 0);
      CHECK(h_1s.coefficients.size() == 3);
      CHECK(h_1s.exponents.size() == 3);

      CHECK(h_1s.coefficients[0] == 0.15432897);
      CHECK(h_1s.coefficients[1] == 0.53532814);
      CHECK(h_1s.coefficients[2] == 0.44463454);
      CHECK(h_1s.exponents[0] == 3.42525091);
      CHECK(h_1s.exponents[1] == 0.62391373);
      CHECK(h_1s.exponents[2] == 0.16885540);
    }  // H

    SECTION("Check magnesium") {
      const std::vector<gaussian::Shell>& sh_mg = b.shells_for_atom(12);
      CHECK(sh_mg.size() == 5);

      const gaussian::Shell& mg_1s = sh_mg[0];
      CHECK(mg_1s.l == 0);
      CHECK(mg_1s.coefficients.size() == 3);
      CHECK(mg_1s.exponents.size() == 3);

      CHECK(mg_1s.coefficients[0] == 0.1543289673);
      CHECK(mg_1s.coefficients[1] == 0.5353281423);
      CHECK(mg_1s.coefficients[2] == 0.4446345422);
      CHECK(mg_1s.exponents[0] == 299.2374000);
      CHECK(mg_1s.exponents[1] == 54.5064700);
      CHECK(mg_1s.exponents[2] == 14.7515800);

      //

      const gaussian::Shell& mg_2s = sh_mg[1];
      CHECK(mg_2s.l == 0);
      CHECK(mg_2s.coefficients.size() == 3);
      CHECK(mg_2s.exponents.size() == 3);

      CHECK(mg_2s.coefficients[0] == -0.09996722919);
      CHECK(mg_2s.coefficients[1] == 0.39951282610);
      CHECK(mg_2s.coefficients[2] == 0.70011546890);
      CHECK(mg_2s.exponents[0] == 15.1218200);
      CHECK(mg_2s.exponents[1] == 3.5139870);
      CHECK(mg_2s.exponents[2] == 1.1428570);

      //

      const gaussian::Shell& mg_2p = sh_mg[2];
      CHECK(mg_2p.l == 1);
      CHECK(mg_2p.coefficients.size() == 3);
      CHECK(mg_2p.exponents.size() == 3);

      CHECK(mg_2p.coefficients[0] == 0.1559162750);
      CHECK(mg_2p.coefficients[1] == 0.6076837186);
      CHECK(mg_2p.coefficients[2] == 0.3919573931);
      CHECK(mg_2p.exponents[0] == 15.1218200);
      CHECK(mg_2p.exponents[1] == 3.5139870);
      CHECK(mg_2p.exponents[2] == 1.1428570);

      //

      const gaussian::Shell& mg_3s = sh_mg[3];
      CHECK(mg_3s.l == 0);
      CHECK(mg_3s.coefficients.size() == 3);
      CHECK(mg_3s.exponents.size() == 3);

      CHECK(mg_3s.coefficients[0] == -0.2196203690);
      CHECK(mg_3s.coefficients[1] == 0.2255954336);
      CHECK(mg_3s.coefficients[2] == 0.9003984260);
      CHECK(mg_3s.exponents[0] == 1.3954480);
      CHECK(mg_3s.exponents[1] == 0.3893260);
      CHECK(mg_3s.exponents[2] == 0.1523800);

      //

      const gaussian::Shell& mg_3p = sh_mg[4];
      CHECK(mg_3p.l == 1);
      CHECK(mg_3p.coefficients.size() == 3);
      CHECK(mg_3p.exponents.size() == 3);

      CHECK(mg_3p.coefficients[0] == 0.01058760429);
      CHECK(mg_3p.coefficients[1] == 0.59516700530);
      CHECK(mg_3p.coefficients[2] == 0.46200101200);
      CHECK(mg_3p.exponents[0] == 1.3954480);
      CHECK(mg_3p.exponents[1] == 0.3893260);
      CHECK(mg_3p.exponents[2] == 0.1523800);
    }  // Mg

    // TODO Check pure vs cartesian for d, f and further!
  }  // sto-3g

  SECTION("Failing Gaussan94 examples") {
    using namespace gaussian;

    BasisSetFileFormat g94 = BasisSetFileFormat::Gaussian94;
    std::stringstream ss("");

    SECTION("Empty string") {
      CHECK_THROWS_AS(read_basisset(ss, g94), ExcInvalidBasisSetFile);
    }  // empty string

    SECTION("Missing ****") {
      std::stringstream ss;
      ss << "!Comment1" << '\n' << "Comment2" << '\n';
      CHECK_THROWS_AS(read_basisset(ss, g94), ExcInvalidBasisSetFile);
    }

    SECTION("Unknown element symbol") {
      std::stringstream ss;
      ss << "****" << '\n'
         << "XX    0" << '\n'
         << "S 1 1.00" << '\n'
         << "1.0 1.0" << '\n'
         << "****" << '\n';
      CHECK_THROWS_AS(read_basisset(ss, g94), ExcInvalidBasisSetFile);
    }

    SECTION("Unknown angular momentum") {
      std::stringstream ss;
      ss << "****" << '\n'
         << "H    0" << '\n'
         << "X 1 1.00" << '\n'
         << "1.0 1.0" << '\n'
         << "****" << '\n';
      CHECK_THROWS_AS(read_basisset(ss, g94), ExcInvalidBasisSetFile);
    }

    SECTION("Unexpected number of primitives") {
      std::stringstream ss;
      ss << "****" << '\n'
         << "H    0" << '\n'
         << "S 2 1.00" << '\n'
         << "1.0 1.0" << '\n'
         << "****" << '\n';
      CHECK_THROWS_AS(read_basisset(ss, g94), ExcInvalidBasisSetFile);
    }

    SECTION("Got S, expected SP") {
      std::stringstream ss;
      ss << "****" << '\n'
         << "H    0" << '\n'
         << "SP 1 1.00" << '\n'
         << "1.0 1.0" << '\n'
         << "****" << '\n';
      CHECK_THROWS_AS(read_basisset(ss, g94), ExcInvalidBasisSetFile);
    }

  }  // failing
}  // read_basisset
}  // namespace tests
}  // namespace gint
