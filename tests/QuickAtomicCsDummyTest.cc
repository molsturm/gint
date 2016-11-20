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
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  const OrbitalType otype = COMPLEX_ATOMIC;
  typedef IntegralLookup<stored_matrix_type, otype> int_lookup_type;

  // The reference data for atomic coulomb sturmians
  // with parameters k = 1, Z = 4, n_max =  3, l_max = 2
  typedef SturmianTestData<stored_matrix_type> refdata_type;

  // Setup parameters for the integral library
  const krims::ParameterMap params{
        {"k_exponent", refdata_type::k_exp},
        {"Z_charge", refdata_type::Z},
        {"n_max", 3},
        {"l_max", 2},
        {"basis_type", "atomic/cs_dummy"},
  };

  IntegralDummyTests<int_lookup_type, refdata_type>::run_all(
        int_lookup_type(params));
}

}  // namespace tests
}  // namespace gint
