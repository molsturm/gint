#pragma once

#include "gint/Integral.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/IntegralCollectionBase.hh"
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

class IntegralCollection: public IntegralCollectionBase<COMPLEX_ATOMIC>{
public:
  typedef IntegralCollectionBase<COMPLEX_ATOMIC> base_type;

  static std::string id, name;

  real_type k_exponent, Z_charge;
  int n_max, l_max;
  sturmint::atomic::cs_dummy::Atomic integral_calculator;

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

  const real_type k, Z;
  const size_t mmax, lmax, nmax;

  // Compute alpha*A*x + beta*y into y
  void apply(const const_multivector_type& x,
             multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1,
             const scalar_type beta  = 0) const override {
    // TODO: This will break due to (m,l,n)-order vs. (n,l,m)-order.
    const real_type *x_ptr = x.data().memptr();
    real_type       *y_ptr = const_cast<real_type*>(y.data().memptr());
    cs::apply_to_full_vectors::nuclear_attraction<real_type>(x_ptr,y_ptr,alpha,beta);
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

  NuclearAttractionIntegralCore(const sturmint::atomic::cs_dummy::Atomic& integral_calculator, real_type k, real_type Z)
    : k(k), Z(Z), 
      mmax(integral_calculator.basis_mmax),
      lmax(integral_calculator.basis_lmax),
      nmax(integral_calculator.basis_nmax),
      m_integral_calculator(integral_calculator){ }

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
  const sturmint::atomic::cs_dummy::Atomic& m_integral_calculator;
};

class OverlapIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef real_stored_mtx_type stored_matrix_type;

  /** \brief Multiplication with a stored matrix */
  // TODO: Change basis order from n,l,m to m,l,n to make multiplication
  // contiguous.
  const size_t mmax, lmax, nmax;
  
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None, const scalar_type c_A = 1,
             const scalar_type c_y = 0) const override {
    const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
    real_type* y_ptr = const_cast<real_type*>(y.data().memptr());

    sturmint::atomic::cs::apply_to_full_vectors::overlap(x_ptr, y_ptr, c_A, c_y,mmax,lmax,nmax);
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
    : mmax(integral_calculator.basis_mmax),
      lmax(integral_calculator.basis_lmax),
      nmax(integral_calculator.basis_nmax),
    m_integral_calculator(integral_calculator) {}

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
  const size_t mmax, lmax, nmax;

  /** \brief Multiplication with a stored matrix */
  void apply(const const_multivector_type& x,
             multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1,
             const scalar_type c_y  = 0) const override {
    
    const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
    real_type* y_ptr = const_cast<real_type*>(y.data().memptr());
    
    sturmint::atomic::cs::apply_to_full_vectors::kinetic(x_ptr,y_ptr,k*k*c_A,c_y);
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
	  mmax(integral_calculator.basis_mmax),
	  lmax(integral_calculator.basis_lmax),
	  nmax(integral_calculator.basis_nmax),
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
  typedef linalgwrap::MultiVector<vector_type>  coefficients_type;
  typedef std::shared_ptr<const coefficients_type> coefficients_ptr_type;  

  bool exchange;  // Is this exchange or Coulomb operator?
  real_type k;    // Exponent scale

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  const int mmax, lmax, nmax;  
  
  /** \brief Multiplication with a stored matrix */
  // J_{aq} = J_{ab} X_{bq} = J_{abcd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{abcd} X_{bq} D_{cd}
  // K_{aq} = K_{ab} X_{bq} = J_{acbd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{acbd} X_{bq} D_{cd}
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1, const scalar_type beta = 0) const override {
    using namespace sturmint::orbital_index;
    using namespace sturmint::atomic::cs_dummy;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    for(size_t i=0;i<y.n_rows();i++)
      for(size_t j=0;j<y.n_cols();j++)
	y(i,j) = (beta != 0? beta*y(i,j) : 0);

    for (size_t a = 0; a < norb; a++)
      for (size_t b = 0; b < norb; b++) {
        real_type JKab = (*this)(a, b);

        for (size_t q = 0; q < x.n_cols(); q++) {
          y(a,q) += alpha * JKab * x(b,q);
        }
      }

    //    const int l_max = m_integral_calculator.lmax1;
    //	  const int n_max = m_integral_calculator.nmax1;
    // //    cout << "\nCalculating "<<(exchange?"exchange":"coulomb")<<" integral
    // //    application.\n";
    // for (size_t b1 = 0; b1 < x.n_elem(); b1++) {
    //   // TODO: Mål, om det bliver hurtigere / nemmere at OpenMP'e,
    //   // af at iterere over b1q samlet.
    //   nlm_t nlm1 = nlmbasis::quantum_numbers_from_index(b1);
    //   int8_t /*n1 = nlm1.n, */ l1 = nlm1.l, m1 = nlm1.m;

    //   for (size_t q = 0; q < x.n_vectors(); q++) {
    //     real_type sum = 0;
    //     real_type y_term = (beta != 0 ? beta * y[q][b1] : 0);

    //     for (size_t b2 = 0; b2 < x.n_elem(); b2++) {
    //       nlm_t nlm2 = nlmbasis::quantum_numbers_from_index(b2);
    //       int8_t /*n2 = nlm2.n, */ l2 = nlm2.l, m2 = nlm2.m;

    //       for (size_t b3 = 0; b3 < n_rows(); b3++) {
    //         nlm_t nlm3 = nlmbasis::quantum_numbers_from_index(b3);
    //         int8_t /*n3 = nlm3.n, l3 = nlm3.l,*/ m3 = nlm3.m;

    //         // TODO: Exchange? Da!
    //         int8_t m4 = m3 - m2 + m1;
    //         int l_parity = (l1 + l2) & 1, m_parity = (m1 + m2) & 1;
    //         int l_min = max(l_parity, ::abs(int(m4)) + ((m_parity + l_parity) & 1));  // TODO: Check.

    //         int B2 = exchange ? b3 : b2;  // Swap b2 and b3 if computing exchange
    //         int B3 = exchange ? b2 : b3;
    //         real_type X_bq = x[q][b2];

    //         for (int8_t l4 = l_min; l4 <= l_max; l4 += 2) {

    //           for (int8_t n4 = l4 + 1; n4 <= n_max; n4++) {
    //             int b4 = nlmbasis::index(nlm_t{n4, l4, m4});

    //             for (size_t p = 0; p < Cocc.n_vectors(); p++)
    //               sum += repulsion[b4 + norb * (B3 + norb * (B2 + norb * b1))] * Cocc[p][b3] * Cocc[p][b4] * X_bq;
    //           }
    //         }
    //       }
    //     }

    //     y[q][b1] = sum + y_term;
    //   }
  }

  scalar_type operator()(size_t a, size_t b) const override {
    return simple1_ver(a,b);
  }
  
  /** \brief return an element of the matrix    */
  // J_{ab} = J_{abcd} Cocc_{cp} Cocc_{dp} = J_{abcd} P_{cd}
  // K_{ab} = J_{cbad} Cocc_{cp} Cocc_{dp} = J_{acbd} P_{cd}

  scalar_type simple2_ver(size_t a, size_t b) const {
    using namespace sturmint::atomic::cs_dummy;
    using namespace sturmint::orbital_index;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    const coefficients_type& Cocc(*coefficients_occupied_ptr);

    nlm_t nlm1 = nlmbasis::quantum_numbers_from_index(a);
    nlm_t nlm2 = nlmbasis::quantum_numbers_from_index(b);
    int8_t /*n1 = nlm1.n,*/ l1 = nlm1.l, m1 = nlm1.m;
    int8_t /*n2 = nlm2.n,*/ l2 = nlm2.l, m2 = nlm2.m;

    real_type sum = 0;

    for (size_t c = 0; c < norb; c++) {
      nlm_t nlm3 = nlmbasis::quantum_numbers_from_index(c);
      int8_t /*n3 = nlm3.n, l3 = nlm3.l,*/ m3 = nlm3.m;

      int8_t m4 = m3 - m2 + m1;
      int l_parity = (l1 + l2) & 1, m_parity = (m1 + m2) & 1;
      int l_min = max(l_parity, ::abs(int(m4)) + ((m_parity + l_parity) & 1));  // TODO: Check!!

      size_t A = exchange ? c : a;  // Swap b and c if computing exchange // TODO: Check
      size_t C = exchange ? a : c;

      assert(a==A && c==C);

      for (int8_t l4 = l_min; l4 <= lmax; l4 += 2) {
        // TODO: Fixed number of n's (same-length dot product for all (l,m)-combinations).
        // Same (l,m), different n's should be contiguous in memory.
        for (int8_t n4 = l4 + 1; n4 <= nmax; n4++) {
          int d = nlmbasis::index(nlm_t{n4, l4, m4});
          real_type density_cd = 0;
          for (size_t p = 0; p < Cocc.n_vectors(); p++)
	    density_cd += Cocc[p][c] * Cocc[p][d];

          sum += repulsion14[d + norb * (C + norb * (b + norb * A))] * density_cd;
        }
      }
    }
    return k*sum;
  }

 
  scalar_type simple1_ver(size_t a, size_t b) const {
    using namespace sturmint::atomic::cs_dummy;
    using namespace sturmint::orbital_index;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    const coefficients_type& Cocc(*coefficients_occupied_ptr);

    real_type sum = 0;

    for (size_t c = 0; c < norb; c++) {
      size_t A = exchange ? c : a;  // Swap b and c if computing exchange // TODO: Check
      size_t C = exchange ? a : c;

      for (size_t d = 0; d < norb; d++) {
        real_type density_cd = 0;
        for (size_t p = 0; p < Cocc.n_vectors(); p++)
	  density_cd += Cocc[p][c] * Cocc[p][d];

	sum += repulsion14[d + norb * (C + norb * (b + norb * A))] * density_cd;
      }
    }
    return k*sum;
    }

    scalar_type static14_ver(size_t a, size_t b) const {
      using namespace sturmint::atomic::cs_dummy;
      assert_greater(a, n_rows());
      assert_greater(b, n_cols());

      // This routine computes
      //
      // J_{ab} = \sum_{cd} P_{cd} < ac | bd >
      //        = \sum_{cd} \sum_{k \in occ} (C^{k})_c (C^{k})_d < ac | bd >
      // where a and b are the same centre, so are c and d
      //
      // or alternatively
      //
      // K_{ab} = \sum_{cd} P_{cd} < ab | cd >
      // where a and c are the same centre, so are b and d

      // Repulsion integrals indexed in shell-pairs
      const auto& i_bbbb = repulsion14;
      const size_t nbas = norb;

      assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

      // Density matrix expression
      auto density_bb = outer_prod_sum(*coefficients_occupied_ptr, *coefficients_occupied_ptr);
      assert_dbg(density_bb.n_rows() == nbas, krims::ExcInternalError());
      assert_dbg(density_bb.n_cols() == nbas, krims::ExcInternalError());

      // Shell pair index for basis functions a and b:
      const size_t ab_pair = a * nbas + b;

      // Sum accumulator variable for this exchange or
      // coulomb  matrix element
      scalar_type mat_ab{0};

      // Double loop over basis functions c and d:
      for (size_t c = 0; c < nbas; ++c) {
        // Shell pair index for basis functions a and c:
        const size_t ac_pair = a * nbas + c;

        for (size_t d = 0; d < nbas; ++d) {
          // Shell pair index for basis functions c and d:
          // or b and d:
          const size_t cd_pair = c * nbas + d;
          const size_t db_pair = d * nbas + b;

          // Perform contraction:
          const scalar_type i_elem = exchange ? i_bbbb[ac_pair * nbas * nbas + db_pair]
                                              : i_bbbb[ab_pair * nbas * nbas + cd_pair];
          mat_ab += density_bb(c, d) * k * i_elem;
        }  // d
      }    // c

      return mat_ab;
    }

    ERICore(const Atomic& integral_calculator, bool exchange, real_type k)
          : exchange(exchange), k(k),
	    mmax(integral_calculator.basis_mmax),
	    lmax(integral_calculator.basis_lmax),
	    nmax(integral_calculator.basis_nmax),
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
      return std::string("Electron Repulsion Integrals, ") + (exchange ? "Exchange" : "Coulomb") +
             " operator";
    }

  private:
    const Atomic& m_integral_calculator;
};

  
  
}  // namespace cs_dummy
}  // namespace atomic
}  // namespace gint
