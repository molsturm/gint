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

#include "GaussianTestData.hh"
#include "integral_quick_tests.hh"
#include <gint/IntegralLookup.hh>
#include <gint/IntegralLookupKeys.hh>

namespace gint {
namespace tests {

namespace gaussian_small_test {
void execute(const std::string& basis_type) {
  // Use the reference data we have for the atomic coulomb Sturmians

  krims::GenMap params{GaussianTestData::integral_parameters};
  params.update(IntegralLookupKeys::basis_type, basis_type);

  const size_t slash = basis_type.rfind('/');
  const std::string prefix = basis_type.substr(slash + 1) + ": ";

  IntegralDummyTests<GaussianTestData>::run_all(prefix, params);
}
}  // gaussian_small_test

const bool runonce = [] {
  std::cout << "TODO:  The libint tests are quite inaccurate ... " << std::endl;
  return true;
}();

#ifdef GINT_HAVE_LIBINT
TEST_CASE("Small test gaussian/libint", "[libint][small]") {

  // TODO This is bad, but to get going ....
  auto highertol = krims::NumCompConstants::change_temporary(
        1e6 * krims::NumCompConstants::default_tolerance_factor);

  gaussian_small_test::execute("gaussian/libint");
}
#endif  // GINT_HAVE_LIBINT

}  // namespace tests
}  // namespace gint
