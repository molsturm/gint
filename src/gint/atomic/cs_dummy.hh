#pragma once

#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/atomic/static14.hh"  // TODO: Included for the horrible hack!
#include "gint/config.hh"

#include <krims/ParameterMap.hh>
#include <sturmint/atomic/cs/cs_atomic.hh>
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <sturmint/atomic/data/cs_dummy.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>

namespace gint {
namespace atomic {
namespace cs_dummy {

#include "gint/real_config.hh"

using namespace sturmint::atomic;
using namespace sturmint::atomic::cs_dummy;

class OverlapIntegralCore;
class NuclearAttractionIntegralCore;
class KineticIntegralCore;
class ERICore;

// This integral class uses (n,l,m)-ordering: {{n,1,nmax},{l,0,n-1},{m,-l,l}}

class IntegralCollection : public IntegralCollectionBase<COMPLEX_ATOMIC> {
 public:
  typedef IntegralCollectionBase<COMPLEX_ATOMIC> base_type;

  static std::string id, name;

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
  IntegralCollection(const krims::ParameterMap& parameters);

  /** Lookup an integral by its identifier string */
  integral_matrix_type lookup_integral(const std::string& integral_id) const override;

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::shared_ptr<base_type> create(const krims::ParameterMap& parameters) {
    return std::make_shared<IntegralCollection>(parameters);
  }
};

// ----------------------------------------------------------------------
//			    INTEGRAL CORES
// ----------------------------------------------------------------------
class NuclearAttractionIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef typename base_type::scalar_type scalar_type;

  const real_type k, Z;
  const size_t nmax;

  bool has_transpose_operation_mode() const override { return true; }

  // Compute alpha*A*x + beta*y into y
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override {
    const real_type* x_ptr = x.data().memptr();
    real_type* y_ptr = const_cast<real_type*>(y.data().memptr());
    cs::apply_to_full_vectors::nuclear_attraction_nlm<real_type>(
          x_ptr, y_ptr, Z * k * c_A, c_y, x.n_cols(), nmax);
  }

  /** \brief return an element of the matrix \f$ {V_0}_{\mu',\mu} = -Zk/n
   * \delta_{\mu',\mu} \f$ */
  scalar_type operator()(size_t row, size_t col) const override {
    using sturmint::orbital_index::nlmbasis;
    if (row != col)
      return 0;
    else {
      int n = nlmbasis::quantum_numbers_from_index(row).n;
      return -Z * k / n;
    }
  }

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
  std::string id() const override { return "atomic/cs_dummy/nuclear_attraction"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Nuclear attraction operator"; }

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
  const size_t nmax;

  bool has_transpose_operation_mode() const override { return true; }
  bool has_apply_inverse() const override { return true; }

  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override {
    const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
    real_type* y_ptr = const_cast<real_type*>(y.data().memptr());

    sturmint::atomic::cs::apply_to_full_vectors::overlap_nlm(x_ptr, y_ptr, c_A, c_y,
                                                             x.n_cols(), nmax);
  }

  void apply_inverse(const const_multivector_type& x, multivector_type& y,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A = 1,
                     const scalar_type c_y = 0) const override {
#ifndef GINT_STATIC_INTEGRALS
    static_assert(false, "Need GINT_STATIC_INTEGRALS=ON for cs_dummy for now.");
#endif
    // TODO: Huge hack, but we don't really want to bother with overlap_inverse apply for
    // nlm-order right
    // now.
    //       How to do that: Compute inverse for each (l,m)-block: jump around in vector.
    using namespace static14;
    const auto& Sinv(detail::Static14Data<stored_mtx_type>::sinv_bb);
    static14::apply_stored_matrix(Sinv, x, y, mode, c_A, c_y);
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override {
    using sturmint::orbital_index::nlmbasis;
    assert_greater(row, n_rows()) assert_greater(col, n_cols());

    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return sturmint::atomic::cs::overlap(mui, muj);
  }

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
  std::string id() const override { return "atomic/cs_dummy/overlap"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Overlap operator"; }

 private:
  const Atomic& m_integral_calculator;
};

class KineticIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef real_stored_mtx_type stored_matrix_type;

  real_type k;  // k-exponent
  const size_t nmax;

  bool has_transpose_operation_mode() const override { return true; }

  /** \brief Multiplication with a stored matrix */
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override {
    const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
    real_type* y_ptr = const_cast<real_type*>(y.data().memptr());

    sturmint::atomic::cs::apply_to_full_vectors::overlap_nlm<scalar_type>(
          x_ptr, y_ptr, (-0.5L * k * k) * c_A, c_y, x.n_cols(), nmax);
    y += c_A * k * k * x;  // kinetic(x) = k^2*x-1/2 overlap(x)
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return k * k * sturmint::atomic::cs::kinetic(mui, muj);
  }

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
  std::string id() const override { return "atomic/cs_dummy/kinetic"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Kinetic energy operator"; }

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
  // J_{aq} = J_{ab} X_{bq} = J_{abcd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{abcd} X_{bq} D_{cd}
  // K_{aq} = K_{ab} X_{bq} = J_{acbd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{acbd} X_{bq} D_{cd}
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1, const scalar_type beta = 0) const override {
    using namespace sturmint::orbital_index;
    using namespace sturmint::atomic::cs_dummy;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    for (size_t i = 0; i < y.n_rows(); i++)
      for (size_t j = 0; j < y.n_cols(); j++) y(i, j) = (beta != 0 ? beta * y(i, j) : 0);

    for (size_t a = 0; a < norb; a++)
      for (size_t b = 0; b < norb; b++) {
        real_type JKab = (*this)(a, b);

        for (size_t q = 0; q < x.n_cols(); q++) {
          y(a, q) += alpha * JKab * x(b, q);
        }
      }
  }

  scalar_type operator()(size_t a, size_t b) const override { return simple1_ver(a, b); }

  /** \brief return an element of the matrix    */
  // J_{ab} = J_{abcd} Cocc_{cp} Cocc_{dp} = J_{abcd} P_{cd}
  // K_{ab} = J_{cbad} Cocc_{cp} Cocc_{dp} = J_{acbd} P_{cd}
  scalar_type simple1_ver(size_t a, size_t b) const {
    using namespace sturmint::atomic::cs_dummy;
    using namespace sturmint::orbital_index;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    const coefficients_type& Cocc(*coefficients_occupied_ptr);

    stored_mtx_type density(norb, norb);
    for (size_t p = 0; p < coefficients_occupied_ptr->n_vectors(); p++) {
      const auto& C = Cocc[p];
      for (size_t c = 0; c < norb; c++)
        for (size_t d = 0; d < norb; d++) density(c, d) += C[c] * C[d];
    }

    real_type sum = 0;
    for (size_t c = 0; c < norb; c++) {
      size_t A = exchange ? c : a;  // Swap b and c if computing exchange // TODO: Check
      size_t C = exchange ? a : c;

      size_t i_abc = norb * (C + norb * (b + norb * A));

      for (size_t d = 0; d < norb; d++) sum += repulsion14[i_abc + d] * density(c, d);
    }
    return k * sum;
  }

  ERICore(const Atomic& integral_calculator, bool exchange, real_type k)
        : exchange(exchange),
          k(k),
          nmax(integral_calculator.nmax),
          m_integral_calculator(integral_calculator) {}

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap                                     */
  virtual void update(const krims::ParameterMap& map) override {
    const std::string occ_coeff_key = Integral<stored_mtx_type>::update_key_coefficients;

    if (!map.exists(occ_coeff_key)) return;

    // Get coefficients as a shared pointer (having ownership)
    coefficients_occupied_ptr = static_cast<coefficients_ptr_type>(
          map.at_ptr<coefficients_type>(occ_coeff_key));

    // We will contract the coefficient row index over the number of
    // basis functions.
    if (coefficients_occupied_ptr->n_vectors() == 0) return;
    assert_size(coefficients_occupied_ptr->n_elem(), m_integral_calculator.n_bas());
  }

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new ERICore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id() const override {
    return std::string("atomic/cs_dummy/ERI_") + (exchange ? "K" : "J");
  }

  /** \brief Get the friendly name of the integral */
  std::string name() const override {
    return std::string("Electron Repulsion Integrals, ") +
           (exchange ? "Exchange" : "Coulomb") + " operator";
  }

 private:
  const Atomic& m_integral_calculator;
};

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace gint
