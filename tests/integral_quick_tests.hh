#pragma once
#include <catch.hpp>
#include <krims/ParameterMap.hh>
#include <linalgwrap/TestingUtils.hh>

namespace gint {
namespace tests {
using namespace linalgwrap;
using namespace krims;

template <typename IntegralLookup, typename RefData>
struct IntegralDummyTests {
  typedef IntegralLookup int_lookup_type;
  typedef typename int_lookup_type::integral_type integral_type;
  typedef typename integral_type::stored_matrix_type stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef RefData data;
  static_assert(
        std::is_same<typename integral_type::stored_matrix_type,
                     typename RefData::stored_matrix_type>::value,
        "Stored matrix types of IntegralLookup and RefData need to agree.");

  static void run_all(const IntegralLookup& integrals) {
    // Obtain integral objects:
    integral_type S_bb = integrals("overlap");
    integral_type T_bb = integrals("kinetic");
    integral_type V0_bb = integrals("nuclear_attraction");
    integral_type J_bb = integrals("coulomb");
    integral_type K_bb = integrals("exchange");

    // The update key we need to update the lazy coulomb and exchange matrices
    const std::string update_key = integral_type::update_key_coefficients;

    // Tolerance levels:
    NumCompAccuracyLevel equality_tol = NumCompAccuracyLevel::Default;
    NumCompAccuracyLevel apply_tol = NumCompAccuracyLevel::Sloppy;
    NumCompAccuracyLevel applyinv_tol = NumCompAccuracyLevel::SuperSloppy;

    SECTION("Test overlap") {
      CHECK(S_bb.is_symmetric());
      REQUIRE(S_bb == numcomp(data::Sref).tolerance(equality_tol));
      CHECK(rc::check("Test apply and extract_block of overlap",
                      make_compare_ref_test(S_bb, data::Sref, apply_tol)));
    }

    SECTION("Test overlap apply_inverse") {
      auto test = [&] {
        auto vec =
              *gen::numeric_tensor<MultiVector<vector_type>>(S_bb.n_rows(), 1)
                     .as("test vector");

        MultiVector<vector_type> tmp1(vec.n_elem(), 1);
        MultiVector<vector_type> tmp2(vec.n_elem(), 1);
        MultiVector<vector_type> res1(vec.n_elem(), 1);
        MultiVector<vector_type> res2(vec.n_elem(), 1);

        S_bb.apply(vec, tmp1);
        S_bb.apply_inverse(vec, tmp2);

        S_bb.apply_inverse(tmp1, res1);
        S_bb.apply(tmp2, res2);

        RC_ASSERT((res1 == numcomp(vec).tolerance(applyinv_tol)));
        RC_ASSERT((res2 == numcomp(vec).tolerance(applyinv_tol)));
      };

      CHECK(rc::check("Test apply_inverse of overlap", test));
    }

    SECTION("Test nuclear attraction") {
      CHECK(V0_bb.is_symmetric());
      REQUIRE(V0_bb == numcomp(data::V0ref).tolerance(equality_tol));
      CHECK(rc::check("Test apply and extract_block of nuclear attraction",
                      make_compare_ref_test(V0_bb, data::V0ref, apply_tol)));
    }

    SECTION("Test kinetic") {
      CHECK(T_bb.is_symmetric());
      REQUIRE(T_bb == numcomp(data::Tref).tolerance(equality_tol));
      CHECK(rc::check("Test apply and extract_block of kinetic",
                      make_compare_ref_test(T_bb, data::Tref, apply_tol)));
    }

    SECTION("Test coulomb: Test case 1") {
      J_bb.update({{update_key,
                    static_cast<coefficients_type>(data::coeffref_bo_1)}});

      CHECK(J_bb.is_symmetric());
      REQUIRE(J_bb == numcomp(data::Jref_for_coeff_1).tolerance(equality_tol));
      CHECK(rc::check(
            "Test apply and extract_block of coulomb 1",
            make_compare_ref_test(J_bb, data::Jref_for_coeff_1, apply_tol)));
    }

    SECTION("Test coulomb: Test case 2") {
      J_bb.update({{update_key,
                    static_cast<coefficients_type>(data::coeffref_bo_2)}});

      CHECK(J_bb.is_symmetric());
      REQUIRE(J_bb == numcomp(data::Jref_for_coeff_2).tolerance(equality_tol));
      CHECK(rc::check(
            "Test apply and extract_block of coulomb 2",
            make_compare_ref_test(J_bb, data::Jref_for_coeff_2, apply_tol)));
    }

    SECTION("Test exchange: Test case 1") {
      K_bb.update({{update_key,
                    static_cast<coefficients_type>(data::coeffref_bo_1)}});

      CHECK(K_bb.is_symmetric());
      REQUIRE(K_bb == numcomp(data::Kref_for_coeff_1).tolerance(equality_tol));
      CHECK(rc::check(
            "Test apply and extract_block of exchange 1",
            make_compare_ref_test(K_bb, data::Kref_for_coeff_1, apply_tol)));
    }

    SECTION("Test exchange: Test case 2") {
      K_bb.update({{update_key,
                    static_cast<coefficients_type>(data::coeffref_bo_2)}});

      CHECK(K_bb.is_symmetric());
      REQUIRE(K_bb == numcomp(data::Kref_for_coeff_2).tolerance(equality_tol));
      CHECK(rc::check(
            "Test apply and extract_block of exchange 2",
            make_compare_ref_test(K_bb, data::Kref_for_coeff_2, apply_tol)));
    }
  }

private:
  // Test generator
  static std::function<void(void)> make_compare_ref_test(
        const integral_type& integral, const stored_matrix_type& ref,
        const NumCompAccuracyLevel tolerance) {
    // TODO use multivector generated by rapidcheck instead of a vector here!
    auto test = [&] {
      auto vec = *gen::numeric_tensor<vector_type>(integral.n_cols());
      RC_ASSERT(((integral * vec) == numcomp(ref * vec).tolerance(tolerance)));

      size_t nrows = *gen::inRange<size_t>(1, ref.n_rows())
                            .as("Number of rows of the extracted matrix");
      size_t ncols = *gen::inRange<size_t>(1, ref.n_cols())
                            .as("Number of cols of the extracted matrix");
      size_t start_row =
            *gen::inRange<size_t>(0, ref.n_rows() - nrows).as("Start row");
      size_t start_col =
            *gen::inRange<size_t>(0, ref.n_cols() - ncols).as("Start col");

      stored_matrix_type xtr_int(nrows, ncols, false);
      stored_matrix_type xtr_ref(nrows, ncols, false);
      integral.extract_block(xtr_int, start_row, start_col);
      ref.extract_block(xtr_ref, start_row, start_col);

      RC_ASSERT((xtr_int == numcomp(xtr_ref).tolerance(tolerance)));
    };
    return test;
  }
};

}  // namespace tests
}  // namespace gint
