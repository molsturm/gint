#include "gint/config.hh"

#include "SturmianTestData.hh"
#include "integral_quick_tests.hh"
#include <gint/IntegralLookup.hh>
#include <linalgwrap/SmallVector.hh>
#include <rapidcheck.h>

namespace gint {
namespace tests {
using namespace linalgwrap;
using namespace krims;

#include "gint/real_config.hh"
  
TEST_CASE("Quick atomic cs_naive test", "[quicktest cs_naive]") {

  typedef typename IntegralLookup<COMPLEX_ATOMIC>::integral_type integral_t;
  //  typedef SturmianTestData<stored_mtx_type> ref_data;
  
  // Setup parameters for the integral library
  const krims::ParameterMap params{
        // {"k_exponent", ref_data::k_exp},
        // {"Z_charge", ref_data::Z},
        {"k_exponent", real_type(1)},
        {"Z_charge",   real_type(2)},
        {"m_max", 2},
        {"l_max", 2},
        {"n_len", 3},
        {"basis_type", "cs_naive"},
  };

  // Compare cs_naive atomic Coulomb Sturmian implementation to
  //         cs_dummy atomic Coulomb Sturmian implementation
  //         (which has been checked against cs_static14).
  IntegralLookup<COMPLEX_ATOMIC> integrals(params);

  // Obtain integral objects:
  integral_t S_bb  = integrals("overlap");
  integral_t T_bb  = integrals("kinetic");
  integral_t V0_bb = integrals("nuclear_attraction");
  integral_t J_bb  = integrals("coulomb");
  integral_t K_bb  = integrals("exchange");
  size_t norb = S_bb.n_cols();
  
  std::ofstream mathfile("/tmp/debug.m");

  auto mout = io::make_formatted_stream_writer<io::Mathematica, double>(mathfile, 1e-12);

  mout.write("S", S_bb);
  mout.write("T", T_bb);
  mout.write("V0",V0_bb);

  // Check apply to identity
  typedef typename stored_mtx_type::vector_type vector_type;
  auto Id = linalgwrap::MultiVector<vector_type>(norb,norb);
  for(size_t i=0;i<norb;i++) Id[i][i] = 1;

  auto SxI  = S_bb*Id;
  auto TxI  = T_bb*Id;
  auto V0xI = V0_bb*Id;

  mout.write("SxI", S_bb);
  mout.write("TxI", T_bb);
  mout.write("V0xI",V0_bb);

  linalgwrap::MultiVector<vector_type> SinvxI(norb,norb);
  S_bb.apply_inverse(Id,SinvxI);

  mout.write("SinvxI",SinvxI);
}

}  // namespace tests
}  // namespace gint
