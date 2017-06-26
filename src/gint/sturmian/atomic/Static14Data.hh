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

#ifdef GINT_HAVE_STATIC_INTEGRALS

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_static14 {
using namespace real_valued;

struct Static14Data {
  /** \name Base integral objects */
  ///@{
  //! Kinetic matrix without the k or Z factors.
  static const stored_matrix_type t_bb_base;

  //! Overlap matrix
  static const stored_matrix_type s_bb;

  //! Inverse of the overlap matrix
  static const stored_matrix_type sinv_bb;

  //! Electron-core interaction matrix without k or Z factors.
  static const stored_matrix_type v0_bb_base;

  //! Two electron integrals without k or Z factors.
  static const stored_matrix_type i_bbbb_base;

  //! All these matrices have exactly 14 rows and cols
  static const size_t nbas = 14;
  ///@}
};

}  // namespace cs_static14
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
#endif  // GINT_HAVE_STATIC_INTEGRALS
