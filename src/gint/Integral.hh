#pragma once
#include "IntegralIdentifier.hh"
#include <krims/SubscriptionPointer.hh>
#include <krims/make_unique.hh>
#include <linalgwrap/LazyMatrix_i.hh>

// TODO: This next header causes the IntegralCoreBase.hh to be part of
//       the public gint interface. Would be good to get rid of that
//       by moving all the implementations of the functions given here
//       to the cc file
#include "IntegralCoreBase.hh"

namespace gint {

// Forward-declare
template <typename StoredMatrix>
class IntegralCoreBase;

template <typename StoredMatrix>
class Integral final : public linalgwrap::LazyMatrix_i<StoredMatrix> {
 public:
  typedef linalgwrap::LazyMatrix_i<StoredMatrix> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;
  typedef IntegralCoreBase<stored_matrix_type> core_type;

  /** \name Constructor, destructor and assignment */
  ///@{
  /** Default constructor: gives rise to an unusable integral object */
  Integral() : m_core_ptr(nullptr) {}

  /** Construct from unique pointer of the integral core.
   *
   * \note That we need a pointer here, since we pass implementations
   * of the core_type */
  Integral(std::unique_ptr<core_type> c) : m_core_ptr{std::move(c)} {
    assert_dbg(m_core_ptr != nullptr, krims::ExcInvalidPointer());
  }

  Integral(const Integral& I) : Integral{I.m_core_ptr->clone()} {
    assert_dbg(m_core_ptr != nullptr, krims::ExcInvalidPointer());
  }

  Integral& operator=(const Integral& I) {
    m_core_ptr = I.m_core_ptr->clone();
    assert_internal(m_core_ptr != nullptr);
    return *this;
  }

  Integral(Integral&&) = default;
  Integral& operator=(Integral&&) = default;
  ~Integral() = default;
  ///@}

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override {
    assert_internal(m_core_ptr != nullptr);
    return m_core_ptr->n_rows();
  }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override {
    assert_internal(m_core_ptr != nullptr);
    return m_core_ptr->n_cols();
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override {
    assert_internal(m_core_ptr != nullptr);
    return m_core_ptr->operator()(row, col);
  }

  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type. **/
  bool has_transpose_operation_mode() const override {
    assert_internal(m_core_ptr != nullptr);
    return m_core_ptr->has_transpose_operation_mode();
  }

  /** Is inverse_apply available for this matrix type */
  bool has_apply_inverse() const override {
    assert_internal(m_core_ptr != nullptr);
    return m_core_ptr->has_apply_inverse();
  }

  /** \brief Compute the Matrix-Multivector application
   *
   * Loosely performs the operation
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   * where mode gives the mode how the current object is to be
   * applied (normal, transposed, conjugate transposed)
   *
   * See the documentation of this function in linalgrwap's
   * LazyMatrixExpression class for details.
   */
  void apply(
        const linalgwrap::MultiVector<
              const linalgwrap::MutableMemoryVector_i<scalar_type>>& x_in,
        linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<scalar_type>>& y_out,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = 1, const scalar_type c_y = 0) const override;

  template <typename VectorIn, typename VectorOut,
            linalgwrap::mat_vec_apply_enabled_t<Integral, VectorIn, VectorOut>...>
  void apply(const linalgwrap::MultiVector<VectorIn>& x,
             linalgwrap::MultiVector<VectorOut>& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = 1, const scalar_type c_y = 0) const {
    using namespace linalgwrap;
    MultiVector<const MutableMemoryVector_i<scalar_type>> x_wrapped(x);
    MultiVector<MutableMemoryVector_i<scalar_type>> y_wrapped(y);
    apply(x_wrapped, y_wrapped, mode, c_this, c_y);
  }

  /** \brief Compute the Inverse-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  template <typename VectorIn, typename VectorOut,
            linalgwrap::mat_vec_apply_enabled_t<Integral, VectorIn, VectorOut>...>
  void apply_inverse(const linalgwrap::MultiVector<VectorIn>& x,
                     linalgwrap::MultiVector<VectorOut>& y,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_this = 1, const scalar_type c_y = 0) const {
    using namespace linalgwrap;
    MultiVector<const MutableMemoryVector_i<scalar_type>> x_wrapped(x);
    MultiVector<MutableMemoryVector_i<scalar_type>> y_wrapped(y);
    apply_inverse(x_wrapped, y_wrapped, mode, c_this, c_y);
  }

  /** \brief Compute the Inverse-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  void apply_inverse(
        const linalgwrap::MultiVector<
              const linalgwrap::MutableMemoryVector_i<scalar_type>>& x_in,
        linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<scalar_type>>& y_out,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = 1, const scalar_type c_y = 0) const override;

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   * where mode gives the mode how the current object is to be
   * applied (normal, transposed, conjugate transposed)
   */
  void mmult(const stored_matrix_type& in, stored_matrix_type& out,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = 1, const scalar_type c_out = 0) const override {
    assert_internal(m_core_ptr != nullptr);
    m_core_ptr->mmult(in, out, mode, c_this, c_out);
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
  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_this = 1,
                     const scalar_type c_M = 0) const override {
    assert_internal(m_core_ptr != nullptr);
    m_core_ptr->extract_block(M, start_row, start_col, mode, c_this, c_M);
  }

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap
   * */
  void update(const krims::GenMap& p) override {
    assert_internal(m_core_ptr != nullptr);
    m_core_ptr->update(p);
  }

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(new Integral<StoredMatrix>(*this));
  }

  /** Identifier of the integral, the structure which is guaranteed
   *  to be unique for each different integral.
   *
   * \note The returned object contains functions to access the type of the
   * integral as human readable strings.
   **/
  IntegralIdentifier id() const {
    assert_internal(m_core_ptr != nullptr);
    return m_core_ptr->id();
  }

 private:
  //! The inner integral core object:
  std::unique_ptr<core_type> m_core_ptr;
};

template <typename IntegralCore, typename... Args>
Integral<typename IntegralCore::stored_matrix_type> make_integral(Args&&... args) {
  typedef typename IntegralCore::stored_matrix_type stored_matrix_type;
  return Integral<stored_matrix_type>(krims::make_unique<IntegralCore>(args...));
}

}  // namespace gint
