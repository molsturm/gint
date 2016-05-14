#pragma once
#include "IntegralCoreBase.hh"
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

  /** \brief Construct a two-electron intregral from
   *  some Integral core type
   */
  Integral(const core_type& c);

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override;

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override;

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& m) const override;

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override;

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap
   * */
  void update(const linalgwrap::ParameterMap& p) override;

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override;

private:
  //! The inner integral core object:
  linalgwrap::SubscriptionPointer<const core_type> m_core_ptr;
};

//
// --------------------------------------------------
//

template <typename StoredMatrix>
Integral<StoredMatrix>::Integral(const core_type& c)
      : m_core_ptr{linalgwrap::make_subscription(c, "Integral")} {
  // TODO The description string should probably be more specific here ...
  assert_dbg(false, linalgwrap::ExcNotImplemented());
}

template <typename StoredMatrix>
typename Integral<StoredMatrix>::size_type Integral<StoredMatrix>::n_rows()
      const {
  assert_dbg(false, linalgwrap::ExcNotImplemented());
  return 0;
}

template <typename StoredMatrix>
typename Integral<StoredMatrix>::size_type Integral<StoredMatrix>::n_cols()
      const {
  assert_dbg(false, linalgwrap::ExcNotImplemented());
  return 0;
}

template <typename StoredMatrix>
typename Integral<StoredMatrix>::stored_matrix_type Integral<StoredMatrix>::
operator*(const stored_matrix_type& m) const {
  assert_dbg(false, linalgwrap::ExcNotImplemented());
  return m;
}

template <typename StoredMatrix>
typename Integral<StoredMatrix>::scalar_type Integral<StoredMatrix>::operator()(
      size_type row, size_type col) const {
  assert_dbg(false, linalgwrap::ExcNotImplemented());
  return 0. * (row - col);
}

template <typename StoredMatrix>
void Integral<StoredMatrix>::update(const linalgwrap::ParameterMap&) {
  assert_dbg(false, linalgwrap::ExcNotImplemented());
}

template <typename StoredMatrix>
typename Integral<StoredMatrix>::lazy_matrix_expression_ptr_type
Integral<StoredMatrix>::clone() const {
  // return a copy enwrapped in the pointer type
  return lazy_matrix_expression_ptr_type(new Integral<StoredMatrix>(*this));
}

}  // namespace gint
