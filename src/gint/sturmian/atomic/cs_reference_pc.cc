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

#include "cs_reference_pc.hh"

#ifdef GINT_HAVE_STURMINT
#include "nlm_order/ERICore.hh"
#include "nlm_order/ERICoreHighPrecision.hh"
#include "nlm_order/ERITensor.hh"
#include "nlm_order/OneElectronIntegralCores.hh"
#include <gint/OrbitalType.hh>
#include <sturmint/atomic/cs/CpxToRealRepulsionCalculator.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_reference_pc {

const std::string IntegralCollection::id = "sturmian/atomic/cs_reference_pc";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : IntegralCollectionBase(parameters) {
  const OrbitalType otype = parameters.at<OrbitalType>(IntegralLookupKeys::orbital_type,
                                                       OrbitalType::REAL_MOLECULAR);

  if (otype == OrbitalType::REAL_ATOMIC) {
    Atomic cpx_calculator{m_system.basis};

    m_integral_calculator_ptr.reset(
          new sturmint::CpxToRealRepulsionCalculator<Atomic>{std::move(cpx_calculator)});
  } else if (otype == OrbitalType::COMPLEX_ATOMIC) {
    m_integral_calculator_ptr.reset(new Atomic(m_system.basis));
  } else {
    assert_throw(false, ExcInvalidIntegralParameters(
                              "OrbitalType may only be REAL_ATOMIC or COMPLEX_ATOMIC"));
  }

  m_eri_tensor_ptr.reset(new ERITensor<Atomic>(*m_integral_calculator_ptr, m_system));
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
      return make_integral<eri_type>(*m_integral_calculator_ptr, m_system,
                                     IntegralIdentifier{id, type});

    default:
      assert_implemented(false);
      return Integral<stored_matrix_type>(nullptr);
  }
}

}  // namespace cs_reference_pc
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
