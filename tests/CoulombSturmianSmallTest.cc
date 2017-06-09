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

#include "SturmianTestData.hh"
#include "integral_quick_tests.hh"
#include <gint/IntegralLookup.hh>
#include <gint/IntegralLookupKeys.hh>

namespace gint {
namespace tests {

namespace cs_small_test {
void execute(const std::string& basis_type) {
  // Use the reference data we have for the atomic coulomb Sturmians

  krims::GenMap params{SturmianTestData::integral_parameters};
  params.update(IntegralLookupKeys::basis_type, basis_type);

  const size_t slash = basis_type.rfind('/');
  const std::string prefix = basis_type.substr(slash + 1) + ": ";

  IntegralDummyTests<SturmianTestData>::run_all(prefix, params);
}
}  // cs_small_test

#ifdef GINT_HAVE_STATIC_INTEGRALS
TEST_CASE("Small test atomic/cs_static14", "[cs_static14][small]") {
  cs_small_test::execute("sturmian/atomic/cs_static14");
}
#endif  // GINT_HAVE_STATIC_INTEGRALS

TEST_CASE("Small test atomic/cs_dummy", "[cs_dummy][small]") {
  cs_small_test::execute("sturmian/atomic/cs_dummy");
}

TEST_CASE("Small test atomic/cs_naive", "[cs_naive][small]") {
  cs_small_test::execute("sturmian/atomic/cs_naive");
}

TEST_CASE("Small test atomic/cs_reference", "[cs_reference][small]") {
  cs_small_test::execute("sturmian/atomic/cs_reference");
}

TEST_CASE("Small test atomic/cs_reference_pc", "[cs_reference_pc][small]") {
  cs_small_test::execute("sturmian/atomic/cs_reference_pc");
}

}  // namespace tests
}  // namespace gint
