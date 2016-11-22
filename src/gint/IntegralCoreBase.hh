#pragma once
#include <krims/ParameterMap.hh>
#include <krims/Subscribable.hh>
#include <linalgwrap/Base/Interfaces.hh>
#include <linalgwrap/MultiVector.hh>

namespace gint {

// TODO we might want to be able to make subscriptions to this thing
template <typename StoredMatrix>
class IntegralCoreBase /* : public linalgwrap::Subscribable */ {
public:
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::size_type size_type;
  typedef typename stored_matrix_type::scalar_type scalar_type;
  typedef typename stored_matrix_type::real_type real_type;

  /** \brief Number of rows of the matrix */
  virtual size_type n_rows() const = 0;

  /** \brief Number of columns of the matrix  */
  virtual size_type n_cols() const = 0;

  /** \brief return an element of the matrix    */
  virtual scalar_type operator()(size_type row, size_type col) const = 0;

  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type
   *
   * These operation modes are important for the functions apply,
   * mmult and extract_block
   **/
  virtual bool has_transpose_operation_mode() const { return false; }

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
        const linalgwrap::MultiVector<
              const linalgwrap::MutableMemoryVector_i<scalar_type>>& x,
        linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<scalar_type>>&
              y,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_y =
              linalgwrap::Constants<scalar_type>::zero) const = 0;

  // TODO Later also provide apply which keeps the type

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
  virtual void mmult(
        const stored_matrix_type& in, stored_matrix_type& out,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_out =
              linalgwrap::Constants<scalar_type>::zero) const {
    // TODO for simplicity we do not force this method to be implemented
    // at the moment
    assert_dbg(false, krims::ExcNotImplemented());
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
  virtual void extract_block(
        stored_matrix_type& M, const size_type start_row,
        const size_type start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M =
              linalgwrap::Constants<scalar_type>::zero) const = 0;

  /** \brief Clone the expression */
  virtual std::unique_ptr<IntegralCoreBase> clone() const = 0;

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap
   * */
  virtual void update(const krims::ParameterMap&) {}

  /** Friendly name of the integral */
  virtual std::string name() const = 0;

  /** Identifier of the integral */
  virtual std::string id() const = 0;
};

}  // namespace gint
