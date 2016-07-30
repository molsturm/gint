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
  stored_matrix_type coeffref_bo{{0.923879532511287, 0.},
                                 {1.306562964876374, 0.},
                                 {0., 0.},
                                 {0., 0.},
                                 {0., 0.919211060789804},
                                 {0.923879532511284, 0.},
                                 {0., 0.},
                                 {0., 0.},
                                 {0., 0.919211060789804},
                                 {0., 0.},
                                 {0., 0.},
                                 {0., 0.},
                                 {0., 0.},
                                 {0., 0.}};

  // Reference Coulomb matrix for those nonlinear coefficients
  stored_matrix_type Jref_for_coeff{
        {1.266286204937442, -0.2568693801546084, 0., 0., 0.,
         -0.1047658029927954, 0., 0., 0., 0., 0., -0.00948398136110219, 0., 0.},
        {-0.2568693801546084, 0.7328637781341365, 0., 0., 0.,
         -0.140793301769866, 0., 0., 0., 0., 0., 0.00875397427897254, 0., 0.},
        {0., 0., 0.844335768660014, 0., 0., 0., -0.07638272034979674, 0., 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0.824217206060768, 0., 0., 0., -0.07764712970593046, 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0.844335768660014, 0., 0., 0., -0.07638272034979674,
         0., 0., 0., 0., 0.},
        {-0.1047658029927954, -0.140793301769866, 0., 0., 0.,
         0.5220627813359466, 0., 0., 0., 0., 0., -0.0006889133959990507, 0.,
         0.},
        {0., 0., -0.07638272034979674, 0., 0., 0., 0.5690079524878392, 0., 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., -0.07764712970593046, 0., 0., 0., 0.5598079730366972, 0.,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0., -0.07638272034979674, 0., 0., 0., 0.5690079524878392,
         0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.6156846969142331, 0., 0., 0.,
         0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.6072961346212196, 0., 0.,
         0.},
        {-0.00948398136110219, 0.00875397427897254, 0., 0., 0.,
         -0.0006889133959990507, 0., 0., 0., 0., 0., 0.6044999471902149, 0.,
         0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.6072961346212196,
         0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0.6156846969142331}};

  // Reference Exchange matrix for those nonlinear coefficients
  stored_matrix_type Kref_for_coeff{
        {0.359888202706249, 0.133732508758423, 0., 0., 0., 0.019648923904622,
         0., 0., 0., 0., 0., -0.0165955980829667, 0., 0.},
        {0.133732508758423, 0.150823395068129, 0., 0., 0., 0.0551389305703789,
         0., 0., 0., 0., 0., 0.00688652633071872, 0., 0.},
        {0., 0., 0.0449844455544027, 0., 0., 0., 0.0164915045496037, 0., 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0.0291566394651632, 0., 0., 0., 0.0109363386834629, 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0., 0.238752182931375, 0., 0., 0., 0.135019585984146, 0.,
         0., 0., 0., 0.},
        {0.019648923904622, 0.0551389305703789, 0., 0., 0., 0.043001713255703,
         0., 0., 0., 0., 0., 0.00491206850843349, 0., 0.},
        {0., 0., 0.0164915045496037, 0., 0., 0., 0.0217152655006732, 0., 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0.0109363386834629, 0., 0., 0., 0.0168060425595382, 0., 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0., 0.135019585984146, 0., 0., 0., 0.103510480185157, 0.,
         0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0114998765395547, 0., 0., 0.,
         0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0089327437574228, 0., 0.,
         0.},
        {-0.0165955980829667, 0.00688652633071872, 0., 0., 0.,
         0.00491206850843349, 0., 0., 0., 0., 0., 0.0104791419166949, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0161390710173709,
         0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0.0259125310594508}};

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

  SECTION("Test exchange") {
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
