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
#include <gint/config.hh>
#ifdef GINT_HAVE_STURMINT

#include <catch.hpp>
#include <gint/IntegralLookup.hh>
#include <gint/IntegralUpdateKeys.hh>
#include <gint/OrbitalType.hh>
#include <gint/Structure.hh>
#include <gint/find_data_file.hh>
#include <krims/DataFiles.hh>
#include <krims/GenMap.hh>
#include <lazyten/MultiVector.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/TestingUtils.hh>
#include <lazyten/io.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>

namespace gint {
namespace tests {
using namespace lazyten;
using namespace krims;
using namespace std;
using namespace sturmint::orbital_index;

void read_test_coefficients(MultiVector<SmallVector<real_type>>& mv, int n_electrons,
                            int nmax, int lmax, int mmax) {
  const std::string datafile = "test-coefficient-" + std::to_string(n_electrons) + "-" +
                               std::to_string(nmax) + "-" + std::to_string(lmax) + "-" +
                               std::to_string(mmax) + ".bin";
  const std::string fullpath = find_data_file(datafile);

  std::vector<real_type> buf;
  read_binary(fullpath, buf, mv.n_elem() * mv.n_vectors());

  for (size_t i = 0, ij = 0; i < mv.n_vectors(); i++) {
    for (size_t j = 0; j < mv.n_elem(); j++, ij++) {
      mv[i][j] = buf[ij];
    }
  }
}

template <typename Integral>
static std::function<void(void)> make_apply_ptr_vector_test(
      const Integral& integral, const typename Integral::stored_matrix_type& ref,
      const NumCompAccuracyLevel tolerance) {
  auto test = [&] {
    typedef typename Integral::stored_matrix_type stored_matrix_type;
    typedef typename stored_matrix_type::scalar_type scalar_type;
    typedef typename stored_matrix_type::vector_type vector_type;
    typedef lazyten::PtrVector<scalar_type> pv_type;

    auto nrow_vecgen = gen::numeric_tensor<vector_type>(integral.n_rows());
    auto vec =
          *gen::with_l2_norm_in_range(0, 2, std::move(nrow_vecgen)).as("test vector");
    auto ptrvec = make_as_multivector<pv_type>(vec.memptr(), vec.size());

    vector_type res_ref(vec.size());
    res_ref(2) = 12.;
    vector_type res(vec.size());
    res(2)           = 12.;
    auto ptr_res_ref = make_as_multivector<pv_type>(res_ref.memptr(), res_ref.size());
    auto ptr_res     = make_as_multivector<pv_type>(res.memptr(), res.size());

    integral.apply(ptrvec, ptr_res, lazyten::Transposed::None, 2., 4.);
    ref.apply(ptrvec, ptr_res_ref, lazyten::Transposed::None, 2., 4.);

    RC_ASSERT((res == numcomp(res_ref).tolerance(tolerance)));
  };
  return test;
}

TEST_CASE("Quick atomic coefficient test", "[quicktest coefficients]") {
  using namespace real_valued;
  typedef IntegralLookup<stored_matrix_type> int_lookup_type;
  typedef Integral<stored_matrix_type> integral_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef const lazyten::MultiVector<const vector_type> coefficients_type;

  int nmax = 5;
  int lmax = 0;
  int mmax = 0;

  unsigned int n_electrons = 2;

  // Setup parameters for the integral library
  const gint::Atom at(n_electrons, {{0, 0, 0}});  // Have a neutral atom
  const krims::GenMap params{{"k_exponent", 1.703},
                             {"structure", gint::Structure{at}},
                             {"n_max", nmax},
                             {"l_max", lmax},
                             {"m_max", mmax},
                             {"orbital_type", OrbitalType::COMPLEX_ATOMIC},
                             {"basis_type", "sturmian/atomic/cs_reference"}};

  vector<nlm_t> nlmbasis(nlmbasis::basis_from_nlm_order(nmax, lmax, mmax));

  int n_orbitals = nlmbasis.size();

  MultiVector<SmallVector<real_type>> Coefs(n_orbitals, n_electrons, false);
  read_test_coefficients(Coefs, n_electrons, nmax, lmax, mmax);

  int_lookup_type integrals(params);

  // Obtain integral objects:
  integral_type S_bb  = integrals.lookup_integral("overlap");
  integral_type T_bb  = integrals.lookup_integral("kinetic");
  integral_type V0_bb = integrals.lookup_integral("nuclear_attraction");
  integral_type J_bb  = integrals.lookup_integral("coulomb");
  integral_type K_bb  = integrals.lookup_integral("exchange");

  // The update key we need to update the lazy coulomb and exchange matrices
  const std::string update_key = IntegralUpdateKeys::coefficients_occupied;

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
    const NumCompAccuracyLevel apply_tol = NumCompAccuracyLevel::Sloppy;

    auto Jbbconv = static_cast<stored_matrix_type>(J_bb);
    CHECK(rc::check("Test J", make_apply_ptr_vector_test(J_bb, Jbbconv, apply_tol)));
  }

  SECTION("Test K") {
    const NumCompAccuracyLevel apply_tol = NumCompAccuracyLevel::Sloppy;

    auto Kbbconv = static_cast<stored_matrix_type>(K_bb);
    CHECK(rc::check("Test K", make_apply_ptr_vector_test(K_bb, Kbbconv, apply_tol)));
  }
}

}  // exchange
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
