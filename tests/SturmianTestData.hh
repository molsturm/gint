#pragma once
#include <gint/config.hh>
#include <krims/GenMap.hh>
#include <krims/Range.hh>
#include <linalgwrap/MultiVector.hh>

namespace gint {
namespace tests {

struct SturmianTestData {
  typedef real_valued::stored_matrix_type stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef real_valued::scalar_type scalar_type;
  typedef linalgwrap::MultiVector<vector_type> mulvector_type;

  //! The parameters to setup the integral library
  static const krims::GenMap integral_parameters;

  static const stored_matrix_type Sref;
  static const stored_matrix_type V0ref;
  static const stored_matrix_type Tref;

  static const mulvector_type coeffref_bo_1;
  static const stored_matrix_type Jref_for_coeff_1;
  static const stored_matrix_type Kref_for_coeff_1;

  static const mulvector_type coeffref_bo_2;
  static const stored_matrix_type Jref_for_coeff_2;
  static const stored_matrix_type Kref_for_coeff_2;

  /** An extracted part of the ERI tensor */
  static const std::vector<scalar_type> Ibbb_extract1;

  /** The corresponding range */
  static const std::array<krims::Range<size_t>, 4> Ibbb_extract1_range;

  static const std::vector<scalar_type> Ibbb_extract2;
  static const std::array<krims::Range<size_t>, 4> Ibbb_extract2_range;

  /** The tensor data obtained if the Electron repulsion tensor
   *  is contacted with the reference coefficient 1 on all places. */
  static const std::vector<scalar_type> Iffff_ref_1111;

  /** The tensor data obtained if the Electron repulsion tensor
   *  is contacted with the reference coefficient 2 on all places. */
  static const std::vector<scalar_type> Iffff_ref_2222;

  /** The tensor data obtained if the Electron repulsion tensor
   *  is contacted with a mixture of the coefficients */
  static const std::vector<scalar_type> Iffff_ref_2121;

  /** The tensor data obtained if the Electron repulsion tensor
   *  is contacted with a mixture of the coefficients */
  static const std::vector<scalar_type> Iffff_ref_2221;
};

}  // namespace tests
}  // namespace gint
