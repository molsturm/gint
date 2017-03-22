#pragma once
#include <catch.hpp>
#include <krims/GenMap.hh>
#include <linalgwrap/TestingUtils.hh>
#include <linalgwrap/io.hh>

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
  typedef typename stored_matrix_type::scalar_type scalar_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef RefData data;
  static_assert(std::is_same<typename integral_type::stored_matrix_type,
                             typename RefData::stored_matrix_type>::value,
                "Stored matrix types of IntegralLookup and RefData need to agree.");

  static void run_all(const std::string& prefix, const IntegralLookup& integrals) {
    // Obtain integral objects:
    integral_type S_bb = integrals.lookup_integral("overlap");
    integral_type T_bb = integrals.lookup_integral("kinetic");
    integral_type V0_bb = integrals.lookup_integral("nuclear_attraction");
    integral_type J_bb = integrals.lookup_integral("coulomb");
    integral_type K_bb = integrals.lookup_integral("exchange");

    // The update key we need to update the lazy coulomb and exchange matrices
    const std::string update_key = integral_type::update_key_coefficients;

    // TODO Tolerance levels (with normed vectors all default
    //    => what happens if we loosen the restriction
    //    => try this!
    NumCompAccuracyLevel symmetric_tol = NumCompAccuracyLevel::Higher;
    NumCompAccuracyLevel equality_tol = NumCompAccuracyLevel::Default;
    NumCompAccuracyLevel apply_tol = NumCompAccuracyLevel::Default;
    NumCompAccuracyLevel applyinv_tol = NumCompAccuracyLevel::Default;

    SECTION(prefix + "Test overlap") {
      // TODO Replace by better function once in krims:
      double hackish_tolerance = numcomp(1.0).tolerance(symmetric_tol).tolerance();
      CHECK(S_bb.is_symmetric(hackish_tolerance));
      REQUIRE((S_bb == numcomp(data::Sref).tolerance(equality_tol)));
      check_apply_to_identity(S_bb, data::Sref, apply_tol);
      CHECK(rc::check(prefix + "Test apply overlap to pointer vectors",
                      make_apply_ptr_vector_test(S_bb, data::Sref, apply_tol)));
      CHECK(rc::check(prefix + "Test apply and extract_block of overlap",
                      make_compare_ref_test(S_bb, data::Sref, apply_tol)));
    }

    SECTION(prefix + "Test overlap apply_inverse") {
      auto test = [&] {
        auto vec = make_as_multivector<vector_type>(
              *gen_normed_vector(S_bb.n_rows()).as("test vector"));

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

      if (S_bb.has_apply_inverse()) {
        CHECK(rc::check(prefix + "Test apply_inverse of overlap", test));
      }
    }

    SECTION(prefix + "Test nuclear attraction") {
      // TODO Replace by better function once in krims:
      double hackish_tolerance = numcomp(1.0).tolerance(symmetric_tol).tolerance();
      CHECK(V0_bb.is_symmetric(hackish_tolerance));
      REQUIRE((V0_bb == numcomp(data::V0ref).tolerance(equality_tol)));
      CHECK(rc::check(prefix + "Test apply nuclear attraction to pointer vectors",
                      make_apply_ptr_vector_test(V0_bb, data::V0ref, apply_tol)));
      CHECK(rc::check(prefix + "Test apply and extract_block of nuclear attraction",
                      make_compare_ref_test(V0_bb, data::V0ref, apply_tol)));
    }

    SECTION(prefix + "Test kinetic") {
      // TODO Replace by better function once in krims:
      double hackish_tolerance = numcomp(1.0).tolerance(symmetric_tol).tolerance();
      CHECK(T_bb.is_symmetric(hackish_tolerance));
      REQUIRE((T_bb == numcomp(data::Tref).tolerance(equality_tol)));
      check_apply_to_identity(T_bb, data::Tref, apply_tol);
      CHECK(rc::check(prefix + "Test apply kinetic to pointer vectors",
                      make_apply_ptr_vector_test(T_bb, data::Tref, apply_tol)));
      CHECK(rc::check(prefix + "Test apply and extract_block of kinetic",
                      make_compare_ref_test(T_bb, data::Tref, apply_tol)));
    }

    SECTION(prefix + "Test coulomb: Test case 1") {
      J_bb.update({{update_key, static_cast<coefficients_type>(data::coeffref_bo_1)}});

      // TODO Replace by better function once in krims:
      double hackish_tolerance = numcomp(1.0).tolerance(symmetric_tol).tolerance();
      CHECK(J_bb.is_symmetric(hackish_tolerance));
      REQUIRE(J_bb == numcomp(data::Jref_for_coeff_1).tolerance(equality_tol));
      REQUIRE((J_bb == numcomp(data::Jref_for_coeff_1).tolerance(equality_tol)));
      check_apply_to_identity(J_bb, data::Jref_for_coeff_1, apply_tol);
      CHECK(rc::check(
            prefix + "Test apply coulomb 1 to pointer vectors",
            make_apply_ptr_vector_test(J_bb, data::Jref_for_coeff_1, apply_tol)));
      CHECK(rc::check(prefix + "Test apply and extract_block of coulomb 1",
                      make_compare_ref_test(J_bb, data::Jref_for_coeff_1, apply_tol)));
    }

    SECTION(prefix + "Test coulomb: Test case 2") {
      J_bb.update({{update_key, static_cast<coefficients_type>(data::coeffref_bo_2)}});

      // TODO Replace by better function once in krims:
      double hackish_tolerance = numcomp(1.0).tolerance(symmetric_tol).tolerance();
      CHECK(J_bb.is_symmetric(hackish_tolerance));
      REQUIRE(J_bb == numcomp(data::Jref_for_coeff_2).tolerance(equality_tol));
      REQUIRE((J_bb == numcomp(data::Jref_for_coeff_2).tolerance(equality_tol)));
      check_apply_to_identity(J_bb, data::Jref_for_coeff_2, apply_tol);
      CHECK(rc::check(
            prefix + "Test apply coulomb 2 to pointer vectors",
            make_apply_ptr_vector_test(J_bb, data::Jref_for_coeff_2, apply_tol)));
      CHECK(rc::check(prefix + "Test apply and extract_block of coulomb 2",
                      make_compare_ref_test(J_bb, data::Jref_for_coeff_2, apply_tol)));
    }

    SECTION(prefix + "Test exchange: Test case 1") {
      K_bb.update({{update_key, static_cast<coefficients_type>(data::coeffref_bo_1)}});

      // TODO Replace by better function once in krims:
      double hackish_tolerance = numcomp(1.0).tolerance(symmetric_tol).tolerance();
      CHECK(K_bb.is_symmetric(hackish_tolerance));
      REQUIRE(K_bb == numcomp(data::Kref_for_coeff_1).tolerance(equality_tol));
      REQUIRE((K_bb == numcomp(data::Kref_for_coeff_1).tolerance(equality_tol)));
      check_apply_to_identity(K_bb, data::Kref_for_coeff_1, apply_tol);
      CHECK(rc::check(
            prefix + "Test apply exchange 1 to pointer vectors",
            make_apply_ptr_vector_test(K_bb, data::Kref_for_coeff_1, apply_tol)));
      CHECK(rc::check(prefix + "Test apply and extract_block of exchange 1",
                      make_compare_ref_test(K_bb, data::Kref_for_coeff_1, apply_tol)));
    }

    SECTION(prefix + "Test exchange: Test case 2") {
      K_bb.update({{update_key, static_cast<coefficients_type>(data::coeffref_bo_2)}});

      // TODO Replace by better function once in krims:
      double hackish_tolerance = numcomp(1.0).tolerance(symmetric_tol).tolerance();
      CHECK(K_bb.is_symmetric(hackish_tolerance));
      REQUIRE(K_bb == numcomp(data::Kref_for_coeff_2).tolerance(equality_tol));
      REQUIRE((K_bb == numcomp(data::Kref_for_coeff_2).tolerance(equality_tol)));
      check_apply_to_identity(K_bb, data::Kref_for_coeff_2, apply_tol);
      CHECK(rc::check(
            prefix + "Test apply exchange 2 to pointer vectors",
            make_apply_ptr_vector_test(K_bb, data::Kref_for_coeff_2, apply_tol)));
      CHECK(rc::check(prefix + "Test apply and extract_block of exchange 2",
                      make_compare_ref_test(K_bb, data::Kref_for_coeff_2, apply_tol)));
    }
  }

 private:
  static rc::Gen<vector_type> gen_normed_vector(size_t n_cols) {
    return gen::map(gen::numeric_tensor<vector_type>(n_cols), [](vector_type&& v) {
      auto nrm = norm_l2(v);
      return (0. == numcomp(nrm).failure_action(NumCompActionType::Return)) ? v : v / nrm;
    });
  }

  // Test generator
  static std::function<void(void)> make_compare_ref_test(
        const integral_type& integral, const stored_matrix_type& ref,
        const NumCompAccuracyLevel tolerance) {
    // TODO use multivector generated by rapidcheck instead of a vector here!
    auto test = [&] {
      auto vec = *gen_normed_vector(integral.n_cols()).as("Input vector");
      auto mvec = make_as_multivector<vector_type>(vec);

      RC_ASSERT(((integral * vec) == numcomp(ref * vec).tolerance(tolerance)));
      RC_ASSERT(((integral * mvec) == numcomp(ref * mvec).tolerance(tolerance)));

      size_t nrows = *gen::inRange<size_t>(1, ref.n_rows())
                            .as("Number of rows of the extracted matrix");
      size_t ncols = *gen::inRange<size_t>(1, ref.n_cols())
                            .as("Number of cols of the extracted matrix");
      size_t start_row = *gen::inRange<size_t>(0, ref.n_rows() - nrows).as("Start row");
      size_t start_col = *gen::inRange<size_t>(0, ref.n_cols() - ncols).as("Start col");

      stored_matrix_type xtr_int(nrows, ncols, false);
      stored_matrix_type xtr_ref(nrows, ncols, false);
      integral.extract_block(xtr_int, start_row, start_col);
      ref.extract_block(xtr_ref, start_row, start_col);

      RC_ASSERT((xtr_int == numcomp(xtr_ref).tolerance(tolerance)));
    };
    return test;
  }

  static std::function<void(void)> make_apply_ptr_vector_test(
        const integral_type& integral, const stored_matrix_type& ref,
        const NumCompAccuracyLevel tolerance) {
    auto test = [&] {
      typedef linalgwrap::PtrVector<scalar_type> pv_type;
      auto vec = *gen_normed_vector(integral.n_rows()).as("test vector");
      auto ptrvec = make_as_multivector<pv_type>(vec.memptr(), vec.size());

      vector_type res_ref(vec.size());
      res_ref(2) = 12.;
      vector_type res(vec.size());
      res(2) = 12.;
      auto ptr_res_ref = make_as_multivector<pv_type>(res_ref.memptr(), res_ref.size());
      auto ptr_res = make_as_multivector<pv_type>(res.memptr(), res.size());

      integral.apply(ptrvec, ptr_res, linalgwrap::Transposed::None, 2., 4.);
      ref.apply(ptrvec, ptr_res_ref, linalgwrap::Transposed::None, 2., 4.);

      RC_ASSERT((res == numcomp(res_ref).tolerance(tolerance)));
    };
    return test;
  }

  static void check_apply_to_identity(const integral_type& integral,
                                      const stored_matrix_type& ref,
                                      const NumCompAccuracyLevel tolerance) {
    INFO("Application of operator to identity:");
    auto Id = linalgwrap::MultiVector<vector_type>(integral.n_cols(), integral.n_cols());
    for (size_t i = 0; i < integral.n_cols(); i++) Id[i][i] = 1;

    auto AxI = integral * Id;
    auto RefxI = ref * Id;

    REQUIRE((AxI == numcomp(RefxI).tolerance(tolerance)));
  }
};

}  // namespace tests
}  // namespace gint
