#include "SturmianTestData.hh"
#include "integral_quick_tests.hh"
#include <gint/IntegralLookup.hh>
#include <linalgwrap/SmallVector.hh>
#include <rapidcheck.h>

namespace gint {
namespace tests {
using namespace linalgwrap;
using namespace krims;

TEST_CASE("Quick atomic cs_dummy test", "[quicktest cs_dummy]") {
  const OrbitalType otype = OrbitalType::COMPLEX_ATOMIC;
  typedef IntegralLookup<otype> int_lookup_type;

  // The reference data for atomic coulomb sturmians
  // with parameters k = 1, Z = 4, n_max =  3, l_max = 2
  typedef SturmianTestData<real_stored_mtx_type> refdata_type;

  // Setup parameters for the integral library
  const krims::GenMap params{
        {"k_exponent", refdata_type::k_exp},
        {"Z_charge", refdata_type::Z},
        {"n_max", 3},
        {"basis_type", "atomic/cs_dummy"},
  };

  IntegralDummyTests<int_lookup_type, refdata_type>::run_all("cs_dummy: ",
                                                             int_lookup_type(params));
}

}  // namespace tests
}  // namespace gint
