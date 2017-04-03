#pragma once
#ifdef GINT_HAVE_STATIC_INTEGRALS
#include "gint/config.hh"

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
