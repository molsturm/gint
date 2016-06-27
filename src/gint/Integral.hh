#pragma once
#include "IntegralCoreBase.hh"
#include "make_unique.hh"
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/SubscriptionPointer.hh>

namespace gint {

template <typename StoredMatrix>
class Integral : public linalgwrap::LazyMatrix_i<StoredMatrix> {
public:
  typedef linalgwrap::LazyMatrix_i<StoredMatrix> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;
  typedef IntegralCoreBase<stored_matrix_type> core_type;

  /** \name Constructor, destructor and assignment */
  ///@{
  /** Construct from unique pointer of the integral core.
   *
   * \note That we need a pointer here, since we pass implementations
   * of the core_type */

  Integral(std::unique_ptr<core_type> c) : m_core_ptr{std::move(c)} {
    assert_dbg(m_core_ptr != nullptr, linalgwrap::ExcInvalidPointer());
  }

  Integral(const Integral& I) : Integral{I.m_core_ptr->clone()} {}

  Integral& operator=(const Integral& I) {
    m_core_ptr = I.m_core_ptr->clone();
    assert_dbg(m_core_ptr != nullptr, linalgwrap::ExcInternalError());
  }

  Integral(Integral&&) = default;
  Integral& operator=(Integral&&) = default;
  ~Integral() = default;
  ///@}

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_core_ptr->n_rows(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_core_ptr->n_cols(); }

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    return m_core_ptr->operator*(X);
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    return m_core_ptr->operator()(row, col);
  }

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap
   * */
  void update(const linalgwrap::ParameterMap& p) override {
    m_core_ptr->update(p);
  }

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(new Integral<StoredMatrix>(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id() const { return m_core_ptr->id(); }

  /** \brief Get the friendly name of the integral */
  std::string name() const { return m_core_ptr->name(); }

private:
  //! The inner integral core object:
  std::unique_ptr<core_type> m_core_ptr;
};

template <typename IntegralCore, typename... Args>
Integral<typename IntegralCore::stored_matrix_type> make_integral(
      Args&&... args) {
  typedef typename IntegralCore::stored_matrix_type stored_matrix_type;
  return Integral<stored_matrix_type>(make_unique<IntegralCore>(args...));
}

}  // namespace gint
