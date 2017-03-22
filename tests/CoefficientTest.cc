#include <catch.hpp>
#include <gint/IntegralLookup.hh>
#include <gint/config.hh>
#include <gint/io.hh>
#include <krims/GenMap.hh>
#include <linalgwrap/MultiVector.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/TestingUtils.hh>
#include <linalgwrap/io.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>

namespace gint {
namespace tests {
using namespace linalgwrap;
using namespace krims;
using namespace std;
using namespace sturmint::orbital_index;

MultiVector<SmallVector<real_type>> multivector_from_pointer(const real_type* data,
                                                             size_t n_rows,
                                                             size_t n_vectors) {
  MultiVector<SmallVector<real_type>> mv(n_rows, n_vectors);

  for (size_t i = 0, ij = 0; i < n_vectors; i++)
    for (size_t j = 0; j < n_rows; j++, ij++) mv[i][j] = data[ij];

  return mv;
}

void multivector_from_vector(MultiVector<SmallVector<real_type>>& mv,
                             const vector<real_type>& data) {
  assert_size(mv.n_elem() * mv.n_vectors(), data.size());

  for (size_t i = 0, ij = 0; i < mv.n_vectors(); i++)
    for (size_t j = 0; j < mv.n_elem(); j++, ij++) mv[i][j] = data[ij];
}

template <typename Vector>
static rc::Gen<Vector> gen_normed_vector(size_t n_cols) {
  return gen::map(gen::numeric_tensor<Vector>(n_cols), [](Vector&& v) {
    auto nrm = norm_l2(v);
    return (0. == numcomp(nrm).failure_action(NumCompActionType::Return)) ? v : v / nrm;
  });
}

template <typename Integral>
static std::function<void(void)> make_apply_ptr_vector_test(
      const Integral& integral, const typename Integral::stored_matrix_type& ref,
      const NumCompAccuracyLevel tolerance) {
  auto test = [&] {
    typedef typename Integral::stored_matrix_type stored_matrix_type;
    typedef typename stored_matrix_type::scalar_type scalar_type;
    typedef typename stored_matrix_type::vector_type vector_type;
    typedef linalgwrap::PtrVector<scalar_type> pv_type;
    auto vec = *gen_normed_vector<vector_type>(integral.n_rows()).as("test vector");
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

TEST_CASE("Quick atomic coefficient test", "[quicktest coefficients]") {
  const OrbitalType otype = OrbitalType::COMPLEX_ATOMIC;
  typedef IntegralLookup<otype> int_lookup_type;
  typedef typename int_lookup_type::integral_type integral_type;
  typedef typename integral_type::stored_matrix_type stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;

  int n_electrons = 2;
  int nmax = 5;
  int lmax = 0;
  int mmax = 0;

  // Setup parameters for the integral library
  const krims::GenMap params{
        {"k_exponent", 1.703}, {"Z_charge", real_type(n_electrons)},  // Neutral atom
        {"n_max", nmax},       {"l_max", lmax},
        {"m_max", mmax},       {"basis_type", "atomic/cs_reference"}};

  vector<nlm_t> nlmbasis(nlmbasis::basis_from_nlm_order(nmax, lmax, mmax));

  int n_orbitals = nlmbasis.size();

  vector<real_type> Coefs_vector(n_electrons * nlmbasis.size());
  read_test_coefficients(Coefs_vector, n_electrons, nmax, lmax, mmax);

  MultiVector<SmallVector<real_type>> Coefs(
        multivector_from_pointer(&Coefs_vector[0], n_orbitals, n_electrons));

  int_lookup_type integrals(params);

  // Obtain integral objects:
  integral_type S_bb = integrals.lookup_integral("overlap");
  integral_type T_bb = integrals.lookup_integral("kinetic");
  integral_type V0_bb = integrals.lookup_integral("nuclear_attraction");
  integral_type J_bb = integrals.lookup_integral("coulomb");
  integral_type K_bb = integrals.lookup_integral("exchange");

  // The update key we need to update the lazy coulomb and exchange matrices
  const std::string update_key = integral_type::update_key_coefficients;

  J_bb.update({{update_key, static_cast<coefficients_type>(Coefs)}});
  K_bb.update({{update_key, static_cast<coefficients_type>(Coefs)}});

  ofstream debug_file("/tmp/CoefficientTest.m");
  auto debug_out = io::make_writer<io::Mathematica, real_type>(debug_file);

  debug_out.write("Sbb", S_bb);
  debug_out.write("V0bb", V0_bb);
  debug_out.write("Tbb", T_bb);
  debug_out.write("Jbb", J_bb);
  debug_out.write("Kbb", K_bb);
  debug_out.write("JbbConv", static_cast<stored_matrix_type>(J_bb));
  debug_out.write("KbbConv", static_cast<stored_matrix_type>(K_bb));
  debug_out.write("Coefs", Coefs);
  debug_out.write("nmax", nmax);
  debug_out.write("lmax", lmax);
  debug_out.write("mmax", mmax);

  SECTION("Test J") {
    auto Jbbconv = static_cast<stored_matrix_type>(J_bb);
    const NumCompAccuracyLevel apply_tol = NumCompAccuracyLevel::Default;
    CHECK(rc::check("Test J", make_apply_ptr_vector_test(J_bb, Jbbconv, apply_tol)));
  }

  SECTION("Test K") {
    auto Kbbconv = static_cast<stored_matrix_type>(K_bb);
    const NumCompAccuracyLevel apply_tol = NumCompAccuracyLevel::Default;
    CHECK(rc::check("Test K", make_apply_ptr_vector_test(K_bb, Kbbconv, apply_tol)));
  }
}

}  // exchange
}  // gint
