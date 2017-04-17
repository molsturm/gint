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
#include <gint/qnum.hh>

namespace gint {
namespace tests {
using namespace gaussian;

template <typename T>
void check_vectors(const std::vector<T>& c, const std::vector<T>& d) {
  CHECK(c.size() == d.size());
  for (size_t i = 0; i < c.size(); ++i) {
    CHECK(c[i] == d[i]);
  }
}

void check_am_pure(const Shell& s, char am_char) {
  CHECK(s.l == qnum::letter_to_am(am_char));
  if (s.l < 2) {
    CHECK_FALSE(s.pure);
  } else {
    CHECK(s.pure);
  }
}

TEST_CASE("read_basisset tests", "[read_basisset]") {
  //
  // Failling
  //
  SECTION("Failing Gaussan94 examples") {
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

  //
  // sto-3g
  //
  SECTION("sto-3g test") {
    std::ifstream f(find_data_file("basis/sto-3g.g94"));
    const BasisSet b = read_basisset(f, BasisSetFileFormat::Gaussian94);

    SECTION("Check hydrogen") {
      const std::vector<Shell>& h = b.shells_for_atom(1);
      CHECK(h.size() == 1);

      check_am_pure(h[0], 's');
      check_vectors(h[0].coefficients, {0.15432897, 0.53532814, 0.44463454});
      check_vectors(h[0].exponents, {3.42525091, 0.62391373, 0.16885540});
    }  // H

    SECTION("Check magnesium") {
      const std::vector<Shell>& mg = b.shells_for_atom(12);
      CHECK(mg.size() == 5);

      check_am_pure(mg[0], 's');
      check_vectors(mg[0].coefficients, {0.1543289673, 0.5353281423, 0.4446345422});
      check_vectors(mg[0].exponents, {299.2374000, 54.5064700, 14.7515800});

      check_am_pure(mg[1], 's');
      check_vectors(mg[1].coefficients, {-0.09996722919, 0.39951282610, 0.70011546890});
      check_vectors(mg[1].exponents, {15.1218200, 3.5139870, 1.1428570});

      check_am_pure(mg[2], 'p');
      check_vectors(mg[2].coefficients, {0.1559162750, 0.6076837186, 0.3919573931});
      check_vectors(mg[2].exponents, {15.1218200, 3.5139870, 1.1428570});

      check_am_pure(mg[3], 's');
      check_vectors(mg[3].coefficients, {-0.2196203690, 0.2255954336, 0.9003984260});
      check_vectors(mg[3].exponents, {1.3954480, 0.3893260, 0.1523800});

      check_am_pure(mg[4], 'p');
      check_vectors(mg[4].coefficients, {0.01058760429, 0.59516700530, 0.46200101200});
      check_vectors(mg[4].exponents, {1.3954480, 0.3893260, 0.1523800});
    }  // Mg
  }    // sto-3g

  //
  // 6-31g*
  //
  SECTION("6-31g* test") {
    std::ifstream f(find_data_file("basis/6-31g*.g94"));
    const BasisSet b = read_basisset(f, BasisSetFileFormat::Gaussian94);

    SECTION("Check hydrogen") {
      const std::vector<Shell>& h = b.shells_for_atom(1);
      CHECK(h.size() == 2);

      INFO("Checking h shell " << 0);
      check_am_pure(h[0], 's');
      check_vectors(h[0].exponents, {18.7311370, 2.8253937, 0.6401217});
      check_vectors(h[0].coefficients, {0.03349460, 0.23472695, 0.81375733});

      INFO("Checking h shell " << 1);
      check_am_pure(h[1], 's');
      check_vectors(h[1].exponents, {0.1612778});
      check_vectors(h[1].coefficients, {1.0000000});
    }  // H

    SECTION("Check magnesium") {
      const std::vector<Shell>& mg = b.shells_for_atom(12);
      CHECK(mg.size() == 8);

      INFO("Checking mg shell " << 0);
      check_am_pure(mg[0], 's');
      check_vectors(mg[0].exponents, {11722.8000000, 1759.9300000, 400.8460000,
                                      112.8070000, 35.9997000, 12.1828000});
      check_vectors(mg[0].coefficients,
                    {0.0019778, 0.0151140, 0.0739110, 0.2491910, 0.4879280, 0.3196620});

      INFO("Checking mg shell " << 1);
      check_am_pure(mg[1], 's');
      check_vectors(mg[1].exponents, {189.1800000, 45.2119000, 14.3563000, 5.1388600,
                                      1.9065200, 0.7058870});
      check_vectors(mg[1].coefficients, {-0.0032372, -0.0410080, -0.1126000, 0.1486330,
                                         0.6164970, 0.3648290});

      INFO("Checking mg shell " << 2);
      check_am_pure(mg[2], 'p');
      check_vectors(mg[2].exponents, {189.1800000, 45.2119000, 14.3563000, 5.1388600,
                                      1.9065200, 0.7058870});
      check_vectors(mg[2].coefficients,
                    {0.0049281, 0.0349890, 0.1407250, 0.3336420, 0.4449400, 0.2692540});

      INFO("Checking mg shell " << 3);
      check_am_pure(mg[3], 's');
      check_vectors(mg[3].exponents, {0.9293400, 0.2690350, 0.1173790});
      check_vectors(mg[3].coefficients, {-0.2122900, -0.1079850, 1.1758400});

      INFO("Checking mg shell " << 4);
      check_am_pure(mg[4], 'p');
      check_vectors(mg[4].exponents, {0.9293400, 0.2690350, 0.1173790});
      check_vectors(mg[4].coefficients, {-0.0224190, 0.1922700, 0.8461810});

      INFO("Checking mg shell " << 5);
      check_am_pure(mg[5], 's');
      check_vectors(mg[5].exponents, {0.0421061});
      check_vectors(mg[5].coefficients, {1.0000000});

      INFO("Checking mg shell " << 6);
      check_am_pure(mg[6], 'p');
      check_vectors(mg[6].exponents, {0.0421061});
      check_vectors(mg[6].coefficients, {1.0000000});

      INFO("Checking mg shell " << 7);
      check_am_pure(mg[7], 'd');
      check_vectors(mg[7].exponents, {0.1750000});
      check_vectors(mg[7].coefficients, {1.0000000});
    }  // Mg
  }    // 6-31g*

}  // read_basisset
}  // namespace tests
}  // namespace gint
