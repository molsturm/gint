#pragma once

#include "gint/Integral.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/IntegralCollectionBase.hh"
#include <krims/ParameterMap.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>

//TODO: Move all sturmint-specific stuff to sturmint
#include <sturmint/common/common.hh>
#include <sturmint/common/simpletensor.hh>
#include <sturmint/harmonic/Y3Coupling.hh>
#include <sturmint/atomic/data/cs_data.hh>
#include <sturmint/atomic/data/Jlnn-table.hh> // TODO: Organization


namespace gint {
namespace atomic {
namespace cs_naive {

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
//TODO: No! This implementation uses mln-order, but cs_dummy uses nlm-order.
  using gint::atomic::cs_dummy::OverlapIntegralCore;
  using gint::atomic::cs_dummy::KineticIntegralCore;
  using gint::atomic::cs_dummy::NuclearAttractionIntegralCore;
  
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

  
  /** \brief Multiplication with a stored matrix */
  // J_{aq} = J_{ab} X_{bq} = J_{abcd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{abcd} X_{bq} D_{cd}
  // K_{aq} = K_{ab} X_{bq} = J_{acbd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{acbd} X_{bq} D_{cd}
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type alpha = 1, const scalar_type beta = 0) const override {
    using namespace sturmint::orbital_index;
    using namespace sturmint::atomic::cs;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    for (auto& vec : y) linalgwrap::detail::scale_or_set(vec, beta);

    for (size_t a = 0; a < norb; a++)
      for (size_t b = 0; b < norb; b++) {
        real_type JKab = (*this)(a, b);

        for (size_t q = 0; q < x.n_vectors(); q++) {
          y[q][a] += alpha * JKab * x[q][b];
        }
      }
  
  /** \brief return an element of the matrix    */
  // J_{ab} = J_{abcd} Cocc_{cp} Cocc_{dp} = J_{abcd} P_{cd}
  // K_{ab} = J_{cbad} Cocc_{cp} Cocc_{dp} = J_{acbd} P_{cd}

 
  scalar_type operator(size_t a, size_t b) const override {
    using namespace sturmint::atomic::cs_dummy;
    using namespace sturmint::orbital_index;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    const coefficients_type& Cocc(*coefficients_occupied_ptr);
    stored_mtx_type density(norb,norb);
    for (size_t p = 0; p < coefficients_occupied_ptr->n_vectors(); p++){
      const auto& C = Cocc[p];
      for(size_t c=0;c<norb;c++)
	for(size_t d=0;d<norb;d++)
	  density(c,d) += C[c]*C[d];
    }
    
    real_type sum = 0;

    for (size_t c = 0; c < norb; c++) {
      size_t A = exchange ? c : a;  // Swap b and c if computing exchange // TODO: Check
      size_t C = exchange ? a : c;

      for (size_t d = 0; d < norb; d++) {
        real_type density_cd = 0;
        for (size_t p = 0; p < Cocc.n_vectors(); p++)
	  density_cd += Cocc[p][c] * Cocc[p][d];

        sum += sturmint::atomic::cs::repulsion(A,b,C,d) * density_cd;
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
      const auto& i_bbbb = repulsion;
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
          : exchange(exchange), k(k), m_integral_calculator(integral_calculator) {}

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

  
  
}  // namespace cs_naive
}  // namespace atomic
}  // namespace gint
