#include "../../linalgwrap/tests/NumComp.hh"  // TODO This is so quick and dirty ...
#include <catch.hpp>
#include <gint/IntegralLookup.hh>
#include <linalgwrap/ParameterMap.hh>
#include <linalgwrap/SmallVector.hh>
#include <rapidcheck.h>

namespace gint {
namespace tests {

TEST_CASE("Quick atomic cs_dummy test", "[quicktest]") {
  typedef double scalar_type;
  typedef linalgwrap::SmallMatrix<scalar_type> stored_matrix_type;
  typedef linalgwrap::SmallVector<scalar_type> vector_type;
  const OrbitalType otype = COMPLEX_ATOMIC;
  typedef gint::IntegralLookup<stored_matrix_type, otype> int_lookup_type;
  typedef typename int_lookup_type::integral_matrix_type int_term_type;

  // Setup parameters for the integral library
  linalgwrap::ParameterMap intparams;
  intparams.update_copy("k_exponent", 1.0);
  intparams.update_copy("Z_charge", 4.0);
  intparams.update_copy("n_max", 3);
  intparams.update_copy("l_max", 2);

  // Setup integral lookup:
  const std::string basis_type("atomic/cs_dummy");
  int_lookup_type integrals(basis_type, intparams);

  // Obtain integral objects:
  int_term_type S_bb = integrals("overlap");
  int_term_type T_bb = integrals("kinetic");
  int_term_type V0_bb = integrals("nuclear_attraction");
  int_term_type J_bb = integrals("coulomb");
  int_term_type K_bb = integrals("exchange");

  // Reference data:
  stored_matrix_type Sref{
        {1, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},             // row1
        {-0.5, 1, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0},          // row2
        {0, 0, 1, 0, 0, 0, -(1 / sqrt(6.)), 0, 0, 0, 0, 0, 0, 0},  // row3
        {0, 0, 0, 1, 0, 0, 0, -(1 / sqrt(6.)), 0, 0, 0, 0, 0, 0},  // row4
        {0, 0, 0, 0, 1, 0, 0, 0, -(1 / sqrt(6.)), 0, 0, 0, 0, 0},  // row5
        {0, -0.5, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},             // row6
        {0, 0, -(1 / sqrt(6.)), 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},  // row7
        {0, 0, 0, -(1 / sqrt(6.)), 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},  // row8
        {0, 0, 0, 0, -(1 / sqrt(6.)), 0, 0, 0, 1, 0, 0, 0, 0, 0},  // row9
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},                // row10
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},                // row11
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},                // row12
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},                // row13
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}                 // row14
  };

  stored_matrix_type V0ref{
        {-4., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},       // 1
        {0., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},       // 2
        {0., 0., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},       // 3
        {0., 0., 0., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},       // 4
        {0., 0., 0., 0., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0.},       // 5
        {0., 0., 0., 0., 0., -4. / 3., 0., 0., 0., 0., 0., 0., 0., 0.},  // 6
        {0., 0., 0., 0., 0., 0., -4. / 3., 0., 0., 0., 0., 0., 0., 0.},  // 7
        {0., 0., 0., 0., 0., 0., 0., -4. / 3., 0., 0., 0., 0., 0., 0.},  // 8
        {0., 0., 0., 0., 0., 0., 0., 0., -4. / 3., 0., 0., 0., 0., 0.},  // 9
        {0., 0., 0., 0., 0., 0., 0., 0., 0., -4. / 3., 0., 0., 0., 0.},  // 10
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4. / 3., 0., 0., 0.},  // 11
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4. / 3., 0., 0.},  // 12
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4. / 3., 0.},  // 13
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4. / 3.}   // 14
  };

  stored_matrix_type Tref{
        {0.5, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},              // row1
        {0.25, 0.5, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0},           // row2
        {0, 0, 0.5, 0, 0, 0, (0.5 / sqrt(6.)), 0, 0, 0, 0, 0, 0, 0},  // row3
        {0, 0, 0, 0.5, 0, 0, 0, (0.5 / sqrt(6.)), 0, 0, 0, 0, 0, 0},  // row4
        {0, 0, 0, 0, 0.5, 0, 0, 0, (0.5 / sqrt(6.)), 0, 0, 0, 0, 0},  // row5
        {0, 0.25, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0},              // row6
        {0, 0, (0.5 / sqrt(6.)), 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0},  // row7
        {0, 0, 0, (0.5 / sqrt(6.)), 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0},  // row8
        {0, 0, 0, 0, (0.5 / sqrt(6.)), 0, 0, 0, 0.5, 0, 0, 0, 0, 0},  // row9
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0},                 // row10
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0},                 // row11
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0},                 // row12
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0},                 // row13
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5}                  // row14
  };

  // Some non-linear coefficients
  stored_matrix_type coeffref_bo{
        {0.92388, 0.}, {1.30656, 0.}, {0., 0.}, {0., 0.},       {0., 0.919211},
        {0.92388, 0.}, {0., 0.},      {0., 0.}, {0., 0.919211}, {0., 0.},
        {0., 0.},      {0., 0.},      {0., 0.}, {0., 0.}};

  // Reference Coulomb matrix for those nonlinear coefficients
  stored_matrix_type Jref_for_coeff{
        {2.532572409874885, -0.5137387603092169, 0., 0., 0.,
         -0.2095316059855909, 0., 0., 0., 0., 0., -0.01896796272220439, 0., 0.},
        {-0.5137387603092169, 1.465727556268273, 0., 0., 0.,
         -0.2815866035397321, 0., 0., 0., 0., 0., 0.01750794855794507, 0., 0.},
        {0., 0., 1.688671537320028, 0., 0., 0., -0.1527654406995935, 0., 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 1.648434412121535, 0., 0., 0., -0.1552942594118609, 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0., 1.688671537320028, 0., 0., 0., -0.1527654406995935, 0.,
         0., 0., 0., 0.},
        {-0.2095316059855909, -0.2815866035397321, 0., 0., 0.,
         1.044125562671893, 0., 0., 0., 0., 0., -0.001377826791998101, 0., 0.},
        {0., 0., -0.1527654406995935, 0., 0., 0., 1.138015904975678, 0., 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., -0.1552942594118609, 0., 0., 0., 1.119615946073395, 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0., -0.1527654406995935, 0., 0., 0., 1.138015904975678, 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 1.231369393828466, 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.214592269242439, 0., 0., 0.},
        {-0.01896796272220439, 0.01750794855794507, 0., 0., 0.,
         -0.001377826791998101, 0., 0., 0., 0., 0., 1.20899989438043, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.214592269242439, 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         1.231369393828466}};

  // Reference Exchange matrix for those nonlinear coefficients
  stored_matrix_type Kref_for_coeff{
        {0.3598882027062493, 0.1337325087584234, 0., 0., 0.,
         0.01964892390462196, 0., 0., 0., 0., 0., -0.01659559808296671, 0., 0.},
        {0.1337325087584234, 0.1508233950681286, 0., 0., 0., 0.0551389305703789,
         0., 0., 0., 0., 0., 0.006886526330718715, 0., 0.},
        {0., 0., 0.04498444555440265, 0., 0., 0., 0.01649150454960366, 0., 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0.02915663946516317, 0., 0., 0., 0.01093633868346285, 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0.2387521829313747, 0., 0., 0., 0.1350195859841461, 0.,
         0., 0., 0., 0.},
        {0.01964892390462197, 0.0551389305703789, 0., 0., 0.,
         0.04300171325570296, 0., 0., 0., 0., 0., 0.004912068508433492, 0., 0.},
        {0., 0., 0.01649150454960366, 0., 0., 0., 0.02171526550067316, 0., 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0.01093633868346285, 0., 0., 0., 0.0168060425595382, 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0.1350195859841461, 0., 0., 0., 0.1035104801851573, 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.01149987653955469, 0., 0., 0.,
         0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0089327437574228, 0., 0.,
         0.},
        {-0.01659559808296671, 0.006886526330718715, 0., 0., 0.,
         0.004912068508433492, 0., 0., 0., 0., 0., 0.01047914191669485, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.01613907101737086,
         0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0.02591253105945082}};

  // Comparator functor
  typedef linalgwrap::tests::NumComp NumComp;
  double thresh = 1e-12;
  const bool verbose_throw = true;

  // Test generator
  auto make_test = [](const int_term_type& integral,
                      const stored_matrix_type& ref, const double& thresh,
                      const bool verbose_throw) {
    // TODO use matrix generated by rapidcheck instead of a vector here!
    auto test = [&] {
      auto vdata = *rc::gen::container<std::vector<scalar_type>>(
                          integral.n_cols(), rc::gen::arbitrary<scalar_type>())
                          .as("Vector data");
      vector_type vec(vdata);

      RC_ASSERT(NumComp::is_equal_matrix((integral * vec), (ref * vec), thresh,
                                         true, verbose_throw));
    };
    return test;
  };

  // Parameter map containing coeff
  linalgwrap::ParameterMap map;
  map.update_copy("coefficients_occupied", coeffref_bo);

  SECTION("Test overlap") {
    CHECK(S_bb.is_symmetric());
    CHECK(NumComp::is_equal_matrix(S_bb, Sref, thresh, true, verbose_throw));
    CHECK(rc::check("Test application of overlap",
                    make_test(S_bb, Sref, thresh, verbose_throw)));
  }

  SECTION("Test nuclear attraction") {
    CHECK(V0_bb.is_symmetric());
    CHECK(NumComp::is_equal_matrix(V0_bb, V0ref, thresh, true, verbose_throw));
    CHECK(rc::check("Test application of nuclear attraction",
                    make_test(V0_bb, V0ref, thresh, verbose_throw)));
  }

  SECTION("Test kinetic") {
    CHECK(T_bb.is_symmetric());
    CHECK(NumComp::is_equal_matrix(T_bb, Tref, thresh, true, verbose_throw));
    CHECK(rc::check("Test application of kinetic",
                    make_test(T_bb, Tref, thresh, verbose_throw)));
  }

  SECTION("Test coulomb") {
    J_bb.update(map);
    CHECK(J_bb.is_symmetric());
    CHECK(NumComp::is_equal_matrix(J_bb, Jref_for_coeff, thresh, true,
                                   verbose_throw));
    CHECK(rc::check("Test application of coulomb",
                    make_test(J_bb, Jref_for_coeff, thresh, verbose_throw)));
  }

  SECTION("Test coulomb") {
    K_bb.update(map);
    CHECK(K_bb.is_symmetric());
    CHECK(NumComp::is_equal_matrix(K_bb, Kref_for_coeff, thresh, true,
                                   verbose_throw));
    CHECK(rc::check("Test application of exchange",
                    make_test(K_bb, Kref_for_coeff, thresh, verbose_throw)));
  }
}

}  // namespace tests
}  // namespace gint
