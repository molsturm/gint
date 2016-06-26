#pragma once
#include <linalgwrap/Subscribable.hh>
#include <linalgwrap/ParameterMap.hh>

namespace gint {

template <typename StoredMatrix>
class IntegralCoreBase : public linalgwrap::Subscribable { // TODO: Is subscribability necessary / desirable
public:
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::size_type   size_type;
  typedef typename stored_matrix_type::scalar_type scalar_type;

  /** \brief Number of rows of the matrix */
  virtual size_type n_rows() const = 0;

  /** \brief Number of columns of the matrix  */
  virtual size_type n_cols() const = 0;

  /** \brief Multiplication with a stored matrix */
  virtual stored_matrix_type operator*(const stored_matrix_type& X) const = 0;

  /** \brief return an element of the matrix    */
  virtual scalar_type operator()(size_type row, size_type col) const = 0;

  /** \brief Clone the expression */
  virtual std::unique_ptr<IntegralCoreBase> clone() const = 0;

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap
   * */
  virtual void update(const linalgwrap::ParameterMap& p) {}
  
  /** Friendly name of the integral */
  virtual std::string name() const = 0;

  /** Identifier of the integral */
  virtual std::string id() const = 0;
};

}  // namespace gint
