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

#include "cs_naive.hh"

#ifdef GINT_HAVE_STURMINT
#include "gint/sturmian/atomic/NlmBasis.hh"
#include "nlm_order/ERICore.hh"
#include "nlm_order/ERICoreHighPrecision.hh"
#include "nlm_order/ERITensor.hh"
#include "nlm_order/OneElectronIntegralCores.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_naive {

const std::string IntegralCollection::id = "sturmian/atomic/cs_naive";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : IntegralCollectionBase{parameters}, m_integral_calculator{m_system.basis} {
  if (parameters.exists(IntegralLookupKeys::nlm_basis)) {
    const auto& nlmbasis = parameters.at<const NlmBasis>("nlm_basis");

    // Construct an NlmBasis in normal, dense form with nlm order
    const int n_max = nlmbasis.n_max();
    const int l_max = nlmbasis.l_max();
    const int m_max = nlmbasis.m_max();
    const NlmBasis nlmbasis_compare(n_max, l_max, m_max);

    // Implement the general case some day. Most importantly think about the required
    // range checks.
    assert_implemented(nlmbasis_compare == nlmbasis);
  }

  // Check range of n,l,m values
  const int n_max = parameters.at<int>("n_max");
  const int l_max = parameters.at<int>("l_max");
  const int m_max = parameters.at<int>("m_max");
  assert_throw(0 < n_max && n_max <= 6,
               ExcInvalidIntegralParameters(
                     "Maximum principle quantum number (" + std::to_string(n_max) +
                     ") needs to be in the range [1,6] for cs_naive, since higher "
                     "values are not yet implemented."));
  assert_throw(0 <= l_max && l_max <= 4,
               ExcInvalidIntegralParameters(
                     "Maximum angular momentum quantum number (" + std::to_string(l_max) +
                     ") needs to be in the range [0,4] for cs_naive, since higher values "
                     "are not yet implemented."));
  assert_throw(
        0 <= m_max && m_max <= 4,
        ExcInvalidIntegralParameters("Maximum magnetic momentum quantum number (" +
                                     std::to_string(m_max) +
                                     ") needs to be in the range [0,4] for cs_naive, "
                                     "since higher values are not yet implemented."));

  m_eri_tensor_ptr.reset(new ERITensor<Atomic>(m_integral_calculator, m_system));
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  typedef ERICoreHighPrecision<Atomic> eri_type;
  const std::string& id = IntegralCollection::id;

  switch (type) {
    case IntegralType::nuclear_attraction:
      return make_integral<NuclearAttractionIntegralCore>(m_system, id);
    case IntegralType::overlap:
      return make_integral<OverlapIntegralCore>(m_system, id);
    case IntegralType::kinetic:
      return make_integral<KineticIntegralCore>(m_system, id);

    case IntegralType::coulomb: /* nobreak */
    case IntegralType::exchange:
      return make_integral<eri_type>(m_integral_calculator, m_system,
                                     IntegralIdentifier{id, type});

    default:
      assert_implemented(false);
      return Integral<stored_matrix_type>(nullptr);
  }
}

}  // namespace cs_naive
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
