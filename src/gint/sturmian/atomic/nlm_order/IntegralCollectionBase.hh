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

#pragma once
#include "gint/config.hh"

#ifdef GINT_HAVE_STURMINT
#include "SturmintSystem.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralLookupKeys.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

struct IntegralLookupKeys : public gint::IntegralLookupKeys {
  static const std::string k_exponent;
  static const std::string nlm_basis;

  /** TODO Add all of them */
};

class IntegralCollectionBase : public gint::IntegralCollectionBase<stored_matrix_type> {
 public:
  typedef gint::IntegralCollectionBase<stored_matrix_type> base_type;

  /** Parse the parameters and setup the SturmintSystem object. */
  IntegralCollectionBase(const krims::GenMap& parameters);

 protected:
  /** The system information in a way usable by sturmint integrals */
  SturmintSystem m_system;
};

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STURMINT
