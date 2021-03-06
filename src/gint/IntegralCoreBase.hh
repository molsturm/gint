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
#include "IntegralIdentifier.hh"
#include "gint/config.hh"
#include <krims/GenMap.hh>
#include <krims/Subscribable.hh>
#include <krims/TypeUtils.hh>
#include <lazyten/Base/Interfaces.hh>

namespace gint {

// TODO Why does this exist? Could one not use a LazyMatrixExpression all along?

template <typename StoredMatrix>
class IntegralCoreBase {
 public:
  using stored_matrix_type = StoredMatrix;
  using scalar_type        = typename stored_matrix_type::scalar_type;

  static constexpr bool isReal = std::is_same<scalar_type, real_type>::value;
  using multivector_type = krims::conditional_t<isReal, real_valued::multivector_type,
                                                complex_valued::multivector_type>;
  using const_multivector_type =
        krims::conditional_t<isReal, real_valued::const_multivector_type,
                             complex_valued::const_multivector_type>;

  IntegralCoreBase()                        = default;
  virtual ~IntegralCoreBase()               = default;
  IntegralCoreBase(const IntegralCoreBase&) = default;
  IntegralCoreBase(IntegralCoreBase&&)      = default;
  IntegralCoreBase& operator=(const IntegralCoreBase&) = default;
  IntegralCoreBase& operator=(IntegralCoreBase&&) = default;

  /** \brief Number of rows of the matrix */
  virtual size_t n_rows() const = 0;

  /** \brief Number of columns of the matrix  */
  virtual size_t n_cols() const = 0;

  /** \brief return an element of the matrix    */
  virtual scalar_type operator()(size_t row, size_t col) const = 0;

  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type
   *
   * These operation modes are important for the functions apply,
   * mmult and extract_block
   **/
  virtual bool has_transpose_operation_mode() const { return false; }

  /** Is inverse_apply available for this matrix type */
  virtual bool has_apply_inverse() const { return false; }

  /** \brief Compute the Matrix-Multivector application
   *
   * Loosely performs the operation
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   * where mode gives the mode how the current object is to be
   * applied (normal, transposed, conjugate transposed)
   *
   * See the documentation of this function in linalgrwap's
   * LazyMatrixExpression class for details.
   *
   * \note Later we might use an interface which is generic in the
   * vector type, but we're not there yet.
   */
  virtual void apply(
        // NB: This will change when the new multivector interface is implemented.
        const const_multivector_type& x, multivector_type& y,
        const lazyten::Transposed mode = lazyten::Transposed::None,
        const scalar_type c_this = 1, const scalar_type c_y = 0) const;

  /** \brief Compute the Inverse-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  virtual void apply_inverse(
        const const_multivector_type& /*x*/, multivector_type& /*y*/,
        const lazyten::Transposed /*mode*/ = lazyten::Transposed::None,
        const scalar_type /*c_this*/ = 1, const scalar_type /*c_y*/ = 0) const {
    assert_throw(false, krims::ExcDisabled("The apply_inverse function is in general "
                                           "very expensive and is only implemented in "
                                           "some cases. Use the function "
                                           "has_apply_inverse() to check when."));
  }

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   * where mode gives the mode how the current object is to be
   * applied (normal, transposed, conjugate transposed)
   *
   * See the documentation of this function in linalgrwap's
   * LazyMatrixExpression class for details.
   */
  virtual void mmult(const stored_matrix_type& in, stored_matrix_type& out,
                     const lazyten::Transposed mode = lazyten::Transposed::None,
                     const scalar_type c_this = 1, const scalar_type c_out = 0) const {
    // TODO for simplicity we do not force this method to be implemented
    // at the moment
    assert_implemented(false);
    (void)in;
    (void)out;
    (void)c_this;
    (void)c_out;
    (void)mode;
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   *
   *  Loosely speaking we perform
   *  \[ M = c_M \cdot M + (A^{mode})_{rowrange,colrange} \]
   *  where
   *    - rowrange = [start_row, start_row+in.n_rows() ) and
   *    - colrange = [start_col, start_col+in.n_cols() )
   *
   * See the documentation of this function in linalgrwap's
   * LazyMatrixExpression class for details.
   */
  virtual void extract_block(stored_matrix_type& M, const size_t start_row,
                             const size_t start_col,
                             const lazyten::Transposed mode = lazyten::Transposed::None,
                             const scalar_type c_A = 1, const scalar_type c_M = 0) const;

  /** \brief Clone the expression */
  virtual std::unique_ptr<IntegralCoreBase> clone() const = 0;

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap
   * */
  virtual void update(const krims::GenMap&) {}

  /** Identifier of the integral, the structure which is guaranteed
   *  to be unique for each different integral. */
  virtual IntegralIdentifier id() const = 0;
};

}  // namespace gint
