#pragma once

#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/atomic/cs_static14.hh"  // TODO: Included for the horrible hack!
#include "gint/config.hh"

#include <sturmint/atomic/cs/cs_atomic.hh>
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <sturmint/atomic/data/cs_dummy.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>

namespace gint {
namespace atomic {
namespace cs_dummy {

using namespace sturmint::atomic;
using namespace sturmint::atomic::cs_dummy;

// In this namespace all things are real:
typedef real_type scalar_type;
typedef real_stored_mtx_type stored_mtx_type;
typedef real_multivector_type multivector_type;
typedef const_real_multivector_type const_multivector_type;

class OverlapIntegralCore;
class NuclearAttractionIntegralCore;
class KineticIntegralCore;
class ERICore;

// This integral class uses (n,l,m)-ordering: {{n,1,nmax},{l,0,n-1},{m,-l,l}}

class IntegralCollection final
      : public IntegralCollectionBase<OrbitalType::COMPLEX_ATOMIC> {
 public:
  typedef IntegralCollectionBase<OrbitalType::COMPLEX_ATOMIC> base_type;

  const static std::string id;

  real_type k_exponent, Z_charge;
  int n_max;
  sturmint::atomic::cs_dummy::Atomic integral_calculator;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - Z_charge (double): The nuclear change of the system
   *   - n_max (int): The maximal principle quantum number
   */
  IntegralCollection(const krims::GenMap& parameters);

  /** Lookup an integral by its type */
  integral_matrix_type lookup_integral(IntegralType type) const override;

  const std::string& basis_id() const override { return id; }
  std::string basis_name() const override {
    return "Dummy implementation of atomic Coulomb Sturmians";
  }

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::unique_ptr<base_type> create(const krims::GenMap& parameters) {
    return krims::make_unique<IntegralCollection>(parameters);
  }
};

//
// Integral Cores
//

class NuclearAttractionIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef typename base_type::scalar_type scalar_type;

  const real_type k, Z;
  const int nmax;

  bool has_transpose_operation_mode() const override { return true; }

  // Compute alpha*A*x + beta*y into y
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** \brief return an element of the matrix \f$ {V_0}_{\mu',\mu} = -Zk/n
   * \delta_{\mu',\mu} \f$ */
  scalar_type operator()(size_t row, size_t col) const override;

  NuclearAttractionIntegralCore(
        const sturmint::atomic::cs_dummy::Atomic& integral_calculator, real_type k,
        real_type Z)
        : k(k),
          Z(Z),
          nmax(integral_calculator.nmax),
          m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new NuclearAttractionIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override {
    return IntegralIdentifier{IntegralCollection::id, IntegralType::nuclear_attraction};
  }

 private:
  const sturmint::atomic::cs_dummy::Atomic& m_integral_calculator;
};

class OverlapIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef real_stored_mtx_type stored_matrix_type;

  /** \brief Multiplication with a stored matrix */
  // TODO: Change basis order from n,l,m to m,l,n to make multiplication
  // contiguous.
  const int nmax;

  bool has_transpose_operation_mode() const override { return true; }
  bool has_apply_inverse() const override { return true; }

  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  void apply_inverse(const const_multivector_type& x, multivector_type& y,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override;

  OverlapIntegralCore(const Atomic& integral_calculator)
        : nmax(integral_calculator.nmax), m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new OverlapIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override {
    return {IntegralCollection::id, IntegralType::overlap};
  }

 private:
  const Atomic& m_integral_calculator;
};

class KineticIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef real_stored_mtx_type stored_matrix_type;

  real_type k;  // k-exponent
  const int nmax;

  bool has_transpose_operation_mode() const override { return true; }

  /** \brief Multiplication with a stored matrix */
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override;

  KineticIntegralCore(const Atomic& integral_calculator, real_type k)
        : k(k),
          nmax(integral_calculator.nmax),
          m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new KineticIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override {
    return {IntegralCollection::id, IntegralType::kinetic};
  }

 private:
  const Atomic& m_integral_calculator;
};

class ERICore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef real_stored_mtx_type stored_matrix_type;
  typedef typename stored_mtx_type::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  bool exchange;  // Is this exchange or Coulomb operator?
  real_type k;    // Exponent scale

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  const int nmax;

  bool has_transpose_operation_mode() const override { return true; }

  /** \brief Multiplication with a stored matrix */
  // J_{aq} = J_{ab} X_{bq} = J_{abcd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{abcd} X_{bq}
  // D_{cd}
  // K_{aq} = K_{ab} X_{bq} = J_{acbd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{acbd} X_{bq}
  // D_{cd}
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1, const scalar_type beta = 0) const override;

  /** \brief return an element of the matrix    */
  // J_{ab} = J_{abcd} Cocc_{cp} Cocc_{dp} = J_{abcd} P_{cd}
  // K_{ab} = J_{cbad} Cocc_{cp} Cocc_{dp} = J_{acbd} P_{cd}
  scalar_type operator()(size_t a, size_t b) const override;

  ERICore(const Atomic& integral_calculator, bool exchange, real_type k)
        : exchange(exchange),
          k(k),
          nmax(integral_calculator.nmax),
          m_integral_calculator(integral_calculator) {}

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap                                     */
  virtual void update(const krims::GenMap& map) override;

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new ERICore(*this));
  }

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override {
    return {IntegralCollection::id,
            (exchange ? IntegralType::exchange : IntegralType::coulomb)};
  }

 private:
  const Atomic& m_integral_calculator;
};

//
// ---------------------------------------------------------------------------
//

inline void NuclearAttractionIntegralCore::apply(const const_multivector_type& x,
                                                 multivector_type& y,
                                                 const linalgwrap::Transposed mode,
                                                 const scalar_type c_A,
                                                 const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  const real_type* x_ptr = x.data().memptr();
  real_type* y_ptr = const_cast<real_type*>(y.data().memptr());
  cs::apply_to_full_vectors::nuclear_attraction_nlm<real_type>(
        x_ptr, y_ptr, Z * k * c_A, c_y, static_cast<int>(x.n_cols()), nmax);
}

inline scalar_type NuclearAttractionIntegralCore::operator()(size_t row,
                                                             size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  using sturmint::orbital_index::nlmbasis;
  if (row != col)
    return 0;
  else {
    const auto urow = static_cast<unsigned short>(row);
    int n = nlmbasis::quantum_numbers_from_index(urow).n;
    return -Z * k / n;
  }
}

// -----

inline void OverlapIntegralCore::apply(const const_multivector_type& x,
                                       multivector_type& y,
                                       const linalgwrap::Transposed mode,
                                       const scalar_type c_A,
                                       const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
  real_type* y_ptr = const_cast<real_type*>(y.data().memptr());

  sturmint::atomic::cs::apply_to_full_vectors::overlap_nlm(
        x_ptr, y_ptr, c_A, c_y, static_cast<int>(x.n_cols()), nmax);
}

inline void OverlapIntegralCore::apply_inverse(const const_multivector_type& x,
                                               multivector_type& y,
                                               const linalgwrap::Transposed mode,
                                               const scalar_type c_A,
                                               const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_rows());
  assert_size(y.n_rows(), n_cols());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
// All modes are same case since we are symmetric and real, so no
// switching over mode.

#ifndef GINT_STATIC_INTEGRALS
  static_assert(false, "Need GINT_STATIC_INTEGRALS=ON for cs_dummy for now.");
#endif
  // TODO: Huge hack, but we don't really want to bother with overlap_inverse apply for
  // nlm-order right
  // now.
  //       How to do that: Compute inverse for each (l,m)-block: jump around in vector.
  using namespace cs_static14;
  const auto& Sinv(detail::Static14Data<stored_mtx_type>::sinv_bb);
  cs_static14::apply_stored_matrix(Sinv, x, y, mode, c_A, c_y);
}

inline scalar_type OverlapIntegralCore::operator()(size_t row, size_t col) const {
  using sturmint::orbital_index::nlmbasis;
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  const unsigned short urow = static_cast<unsigned short>(row);
  const unsigned short ucol = static_cast<unsigned short>(col);
  const nlm_t mui = nlmbasis::quantum_numbers_from_index(urow),
              muj = nlmbasis::quantum_numbers_from_index(ucol);

  return static_cast<scalar_type>(sturmint::atomic::cs::overlap(mui, muj));
}

// -----

inline void KineticIntegralCore::apply(const const_multivector_type& x,
                                       multivector_type& y,
                                       const linalgwrap::Transposed mode,
                                       const scalar_type c_A,
                                       const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
  real_type* y_ptr = const_cast<real_type*>(y.data().memptr());

  const scalar_type c_kkA = static_cast<scalar_type>(-0.5L * k * k * c_A);
  sturmint::atomic::cs::apply_to_full_vectors::overlap_nlm<scalar_type>(
        x_ptr, y_ptr, c_kkA, c_y, static_cast<int>(x.n_cols()), nmax);
  y += c_A * k * k * x;  // kinetic(x) = k^2*x-1/2 overlap(x)
}

inline scalar_type KineticIntegralCore::operator()(size_t row, size_t col) const {
  using sturmint::orbital_index::nlmbasis;
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  const unsigned short urow = static_cast<unsigned short>(row);
  const unsigned short ucol = static_cast<unsigned short>(col);
  const nlm_t mui = nlmbasis::quantum_numbers_from_index(urow),
              muj = nlmbasis::quantum_numbers_from_index(ucol);

  return static_cast<scalar_type>(k * k * sturmint::atomic::cs::kinetic(mui, muj));
}

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace gint
