#pragma once

#include "gint/Integral.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/IntegralCollectionBase.hh"
#include <krims/ParameterMap.hh>
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <sturmint/atomic/data/cs_dummy.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>


namespace gint {
namespace atomic {
namespace cs_dummy {

#include "gint/real_config.hh"
  
using namespace sturmint::atomic::cs_dummy;

class OverlapIntegralCore;
class NuclearAttractionIntegralCore;
class KineticIntegralCore;
class ERICore;

class IntegralCollection: public IntegralCollectionBase<COMPLEX_ATOMIC>{
public:
  typedef IntegralCollectionBase<COMPLEX_ATOMIC> base_type;

  static std::string id, name;

  real_type k_exponent, Z_charge;
  int n_max, l_max;
  Atomic integral_calculator;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - Z_charge (double): The nuclear change of the system
   *   - n_max (int): The maximal principle quantum number
   *   - l_max (int): Maximal angular momentum qn
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
  typedef sturmint::real_t real_type;

  const real_type k, Z;

  // Compute alpha*A*x + beta*y into y
  void apply(const const_multivector_type& x,
             multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1,
             const scalar_type beta  = 0) const override {
    using namespace linalgwrap;

    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;

    // TODO: Add preconditions (dimension check, etc.)
    // TODO: Basis lexicographic order m,l,n is way more efficient than n,l,m.
    // TODO: Make memory layout more friendly to operate on
    //       many vectors.

    for (size_t i = 0; i < x.n_elem(); i++) {
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n;
      
      for (size_t j = 0; j < x.n_vectors(); j++)
	y[j][i] = alpha * (-Z*k/n) * x[j][i] + (beta != 0? beta*y[j][i] : 0);
    }
  }

  /** \brief return an element of the matrix \f$ {V_0}_{\mu',\mu} = -Zk/n
   * \delta_{\mu',\mu} \f$ */
  scalar_type operator()(size_t row, size_t col) const override {
    using sturmint::orbital_index::nlmbasis;
    if (row != col)
      return 0;
    else {
      int n = nlmbasis::quantum_numbers_from_index(row).n;
      return -Z*k/n;
    }
  }

  NuclearAttractionIntegralCore(const Atomic& integral_calculator, real_type k, real_type Z)
        : k(k), Z(Z), m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new NuclearAttractionIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id() const override {
    return "atomic/cs_dummy/nuclear_attraction";
  }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Nuclear attraction operator"; }

private:
  const Atomic& m_integral_calculator;
};

class OverlapIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef real_stored_mtx_type stored_matrix_type;

  /** \brief Multiplication with a stored matrix */
  // TODO: Change basis order from n,l,m to m,l,n to make multiplication
  // contiguous.

  void apply(const const_multivector_type& x,
             multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1,
             const scalar_type beta  = 0) const override {  
    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;

    //        for (auto& vec : y) linalgwrap::detail::scale_or_set(vec, beta);
    // TODO: Add preconditions (dimension check, etc.)
    // TODO: Basis lexicographic order m,l,n is way more efficient than n,l,m.
    // TODO: Make memory layout more friendly to operate on
    //       many vectors.    
    
    const int nmax = m_integral_calculator.nmax1;

    // Apply S^{(l)}_{n',n} = \delta_{n',n} +
    // S^{(l)}_{n,n+1}(\delta_{n',n+1}+\delta_{n'+1,n}).
    // S is block-diagonal w.r.t. l and m, and only terms with |n'-n|<2 are
    // nonzero.
    for (size_t i = 0; i < x.n_elem(); i++) {
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n, l = nlm.l, m = nlm.m;
      size_t i_nminus = nlmbasis::index(nlm_t(n - 1, l, m)),
             i_nplus  = nlmbasis::index(nlm_t(n + 1, l, m));

      // Overlap is symmetric in n,np
      double S_nnl = m_integral_calculator.overlap(n, n + 1, l);

      for (size_t j = 0; j < x.n_vectors(); j++)
        y[j][i] = ((n > l + 1) ? S_nnl * x[j][i_nminus] : 0)
	          + x[j][i] +
                  ((n < nmax) ? S_nnl * x[j][i_nplus] : 0)
	          +
	          (beta != 0 ? beta * y[j][i] : 0);
    }
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override {
    using sturmint::orbital_index::nlmbasis;
    assert_greater(row, n_rows()) assert_greater(col, n_cols());

    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return m_integral_calculator.overlap(mui, muj);
  }

  OverlapIntegralCore(const Atomic& integral_calculator)
        : m_integral_calculator(integral_calculator) {}

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

  /** \brief Multiplication with a stored matrix */
  void apply(const const_multivector_type& x,
             multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1,
             const scalar_type beta  = 0) const override {  
    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;

    const int nmax = m_integral_calculator.nmax1;

    for (size_t i = 0; i < x.n_elem(); i++) {
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n, l = nlm.l, m = nlm.m;

      // T is symmetric in n1,n2 and has the same sparsity pattern as the
      // overlap matrix S.
      size_t i_nminus = nlmbasis::index(nlm_t(n - 1, l, m)),
             i_nplus = nlmbasis::index(nlm_t(n + 1, l, m));
      double T_nnl = -0.5L * m_integral_calculator.overlap(n, n + 1, l);

      for (size_t j = 0; j < x.n_vectors(); j++) {
        y[j][i] = ((n > l + 1) ? T_nnl * x[j][i_nminus] : 0) + 0.5L * x[j][i] +
                  ((n < nmax)  ? T_nnl * x[j][i_nplus]  : 0)
	  	  +
	          (beta != 0 ? beta * y[j][i] : 0);

        y[j][i] *= k * k;
      }
    }
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return k * k * m_integral_calculator.kinetic(mui, muj);
  }

  KineticIntegralCore(const Atomic& integral_calculator, real_type k)
        : k(k), m_integral_calculator(integral_calculator) {}

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
  typedef multivector_type coefficients_type;
  typedef std::shared_ptr<const coefficients_type> coefficients_ptr_type;  

  bool exchange;  // Is this exchange or Coulomb operator?
  real_type k;    // Exponent scale

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  /** \brief Multiplication with a stored matrix */
  // J_{b1,q} = J_{b1,b2} X_{b2,q} = J_{b1,b2,b3,b4} X_{b2,q} Cocc_{b3,p} Cocc_{b4,p} = J_{b1,b2,b3,b4}
  // X_{b2,q} D_{b3,b4}
  // K_{b1,q} = K_{b1,b2} X_{b2,q} = J_{b1,b3,b2,b4} X_{b2,q} Cocc_{b3,p} Cocc_{b4,p} = J_{b1,b3,b2,b4}
  // X_{b2,q} D_{b3,b4}
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1, const scalar_type beta = 0) const override {
    using namespace sturmint::orbital_index;
    using namespace sturmint::atomic::cs_dummy;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    const int l_max = m_integral_calculator.lmax1;
    const int n_max = m_integral_calculator.nmax1;
    const multivector_type& Cocc(*coefficients_occupied_ptr);

    //    cout << "\nCalculating "<<(exchange?"exchange":"coulomb")<<" integral
    //    application.\n";
    for (size_t b1 = 0; b1 < x.n_elem(); b1++) {
      // TODO: Mål, om det bliver hurtigere / nemmere at OpenMP'e,
      // af at iterere over b1q samlet.
      nlm_t nlm1 = nlmbasis::quantum_numbers_from_index(b1);
      int8_t /*n1 = nlm1.n, */ l1 = nlm1.l, m1 = nlm1.m;

      for (size_t q = 0; q < x.n_vectors(); q++) {
        real_type sum = 0;
        real_type y_term = (beta != 0 ? beta * y[q][b1] : 0);

        for (size_t b2 = 0; b2 < x.n_elem(); b2++) {
          nlm_t nlm2 = nlmbasis::quantum_numbers_from_index(b2);
          int8_t /*n2 = nlm2.n, */ l2 = nlm2.l, m2 = nlm2.m;

          for (size_t b3 = 0; b3 < n_rows(); b3++) {
            nlm_t nlm3 = nlmbasis::quantum_numbers_from_index(b3);
            int8_t /*n3 = nlm3.n, l3 = nlm3.l,*/ m3 = nlm3.m;

            // TODO: Exchange? Da!
            int8_t m4 = m3 - m2 + m1;
            int l_parity = (l1 + l2) & 1, m_parity = (m1 + m2) & 1;
            int l_min = max(l_parity, ::abs(int(m4)) + ((m_parity + l_parity) & 1));  // TODO: Check.

            int B2 = exchange ? b3 : b2;  // Swap b2 and b3 if computing exchange
            int B3 = exchange ? b2 : b3;
            real_type X_bq = x[q][b2];

            for (int8_t l4 = l_min; l4 <= l_max; l4 += 2) {

              for (int8_t n4 = l4 + 1; n4 <= n_max; n4++) {
                int b4 = nlmbasis::index(nlm_t{n4, l4, m4});

                for (size_t p = 0; p < Cocc.n_vectors(); p++)
                  sum += repulsion[b4 + norb * (B3 + norb * (B2 + norb * b1))] * Cocc[p][b3] * Cocc[p][b4] * X_bq;
              }
            }
          }
        }

        y[q][b1] = sum + y_term;
      }
    }
  }

  /** \brief return an element of the matrix    */
  // J_{b1,b2} = J_{b1,b2,b3,b4} Cocc_{b3,p} Cocc_{b4,p} = J_{b1,b2,b3,b4} D_{b3,b4}
  // K_{b1,b2} = J_{b1,b3,b2,b4} Cocc_{b3,p} Cocc_{b4,p} = J_{b1,b2,b3,b4} D_{b3,b4}
  scalar_type operator()(size_t b1, size_t b2) const override {
    using namespace sturmint::atomic::cs_dummy;
    using namespace sturmint::orbital_index;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    const int l_max = m_integral_calculator.lmax1;
    const int n_max = m_integral_calculator.nmax1;
    const multivector_type& Cocc(*coefficients_occupied_ptr);

    nlm_t nlm1 = nlmbasis::quantum_numbers_from_index(b1);
    nlm_t nlm2 = nlmbasis::quantum_numbers_from_index(b2);
    int8_t /*n1 = nlm1.n,*/ l1 = nlm1.l, m1 = nlm1.m;
    int8_t /*n2 = nlm2.n, */ l2 = nlm2.l, m2 = nlm2.m;

    real_type sum = 0;

    for (size_t b3 = 0; b3 < Cocc.n_elem(); b3++) {
      nlm_t nlm3 = nlmbasis::quantum_numbers_from_index(b3);
      int8_t /*n3 = nlm3.n, l3 = nlm3.l,*/ m3 = nlm3.m;

      int8_t m4 = m3 - m2 + m1;
      int l_parity = (l1 + l2) & 1, m_parity = (m1 + m2) & 1;
      int l_min = max(l_parity, ::abs(int(m4)) + ((m_parity + l_parity) & 1));  // TODO: Check!!

      int B2 = exchange ? b3 : b2;  // Swap b2 and b3 if computing exchange // TODO: Check
      int B3 = exchange ? b2 : b3;

      for (int8_t l4 = l_min; l4 <= l_max; l4 += 2) {
        // TODO: Fixed number of n's (same-length dot product for all
        // (l,m)-combinations). Same (l,m), different n's should be
        // contiguous in memory.
        for (int8_t n4 = l4 + 1; n4 <= n_max; n4++) {
          int b4 = nlmbasis::index(nlm_t{n4, l4, m4});

          for (size_t p = 0; p < Cocc.n_vectors(); p++)
            sum += repulsion[b4 + norb * (B3 + norb * (B2 + norb * b1))] * Cocc[p][b3] * Cocc[p][b4];
        }
      }
    }
    return sum;
  }

  ERICore(const Atomic& integral_calculator, bool exchange, real_type k)
        : exchange(exchange),
          k(k),
          m_integral_calculator(integral_calculator) {}

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap                                     */
  virtual void update(const krims::ParameterMap& map) override {
    const std::string occ_coeff_key = Integral<stored_mtx_type>::update_key_coefficients;

    if (!map.exists(occ_coeff_key)) return;

    // Get coefficients as a shared pointer (having ownership)
    coefficients_occupied_ptr =
          static_cast<coefficients_ptr_type>(map.at_ptr<coefficients_type>(occ_coeff_key));

    // We will contract the coefficient row index over the number of
    // basis functions.
    if (coefficients_occupied_ptr->n_vectors() == 0) return;
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
