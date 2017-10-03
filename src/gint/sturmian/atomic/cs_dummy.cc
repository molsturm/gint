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

#include "cs_dummy.hh"

#ifdef GINT_HAVE_STURMINT
#include "gint/sturmian/atomic/NlmBasis.hh"
#include "nlm_order/ERICore.hh"
#include "nlm_order/ERICoreHighPrecision.hh"
#include "nlm_order/ERITensor.hh"
#include "nlm_order/OneElectronIntegralCores.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_dummy {

const std::string IntegralCollection::id = "sturmian/atomic/cs_dummy";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : IntegralCollectionBase(parameters) {
  if (parameters.exists("n_max") and not parameters.exists("repulsiondata_filename")) {
    const int nmax = parameters.at<int>("n_max");
    const int lmax = parameters.at<int>("l_max");
    const int mmax = parameters.at<int>("m_max");

    m_repulsiondata_filename = std::string("repulsiondata-nlm-") + to_string(nmax) + "-" +
                               to_string(lmax) + "-" + to_string(mmax) + ".bin";
  } else {
    m_repulsiondata_filename = parameters.at<string>("repulsiondata_filename");
  }
  m_integral_calculator = Atomic(m_system.basis, m_repulsiondata_filename);
  m_eri_tensor_ptr.reset(new ERITensor<Atomic>(m_integral_calculator, m_system));
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  typedef ERICore<Atomic> eri_type;
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

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
