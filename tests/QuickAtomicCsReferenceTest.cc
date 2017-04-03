#include "SturmianTestData.hh"
#include "gint/config.hh"
#include "integral_quick_tests.hh"
#include <gint/IntegralLookup.hh>
#include <gint/OrbitalType.hh>
#include <linalgwrap/SmallVector.hh>
#include <rapidcheck.h>

namespace gint {
namespace tests {
using namespace linalgwrap;
using namespace krims;

TEST_CASE("Quick atomic cs_reference test", "[quicktest cs_reference]") {
  typedef IntegralLookup<real_valued::stored_matrix_type> int_lookup_type;

  // The reference data for atomic coulomb sturmians
  // with parameters k = 1, Z = 4, n_max =  3, l_max = 2
  typedef SturmianTestData<real_valued::stored_matrix_type> refdata_type;

  // Setup parameters for the integral library
  const krims::GenMap params{
        {"k_exponent", refdata_type::k_exp},
        {"Z_charge", refdata_type::Z},
        {"n_max", 3},
        {"l_max", 2},
        {"m_max", 2},
        {"orbital_type", OrbitalType::COMPLEX_ATOMIC},
        {"basis_type", "atomic/cs_reference"},
  };

  IntegralDummyTests<int_lookup_type, refdata_type>::run_all("cs_reference: ",
                                                             int_lookup_type(params));
}

}  // namespace tests
}  // namespace gint
