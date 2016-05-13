#include <catch.hpp>
#include <gint/version.hh>
#include <rapidcheck.h>
#include <sturmint/version.hh>

namespace gint {
namespace tests {

TEST_CASE("Dummy test", "[dummy]") {

  // just to see if linking works
  std::cout << "gint version: " << version::version_string() << std::endl;
  std::cout << "sturmint version: " << sturmint::version::version_string()
            << std::endl;

  auto test = [](int x) { RC_ASSERT(x + x == 2 * x); };
  CHECK(rc::check("Run dummy", test));
}

}  // namespace tests
}  // namespace gint
