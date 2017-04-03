#pragma once
#include "gint/IntegralCollectionBase.hh"
#include "nlm_order.hh"
#include <sturmint/atomic/cs_reference/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_reference {
using namespace nlm_order;
using sturmint::atomic::cs_reference::Atomic;

// This integral class uses (n,l,m)-ordering: {{n,1,nmax},{l,0,n-1},{m,-l,l}}

class IntegralCollection final : public IntegralCollectionBase<stored_matrix_type> {
 public:
  typedef IntegralCollectionBase<stored_matrix_type> base_type;

  const static std::string id;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - Z_charge (double): The nuclear change of the system
   *   - n_max (int): The maximal principle quantum number
   *   - l_max (int): Maximal azimuthal quantum number
   *   - m_max (int): Maximal magnetic quantum number
   *   - nlmbasis (NlmCollection): The precise basis triples (nlm) to use
   */
  IntegralCollection(const krims::GenMap& parameters);

  using base_type::lookup_integral;
  integral_matrix_type lookup_integral(IntegralType type) const override;

  const std::string& basis_id() const override { return id; }
  std::string basis_name() const override {
    return "Reference implementation of atomic Coulomb Sturmians";
  }

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::unique_ptr<base_type> create(const krims::GenMap& parameters) {
    return krims::make_unique<IntegralCollection>(parameters);
  }

 private:
  /** The system information in a way usable by sturmint integrals */
  SturmintSystem m_system;

  /** The integral calculator object */
  Atomic m_integral_calculator;
};

//

class ERICore : public IntegralCoreBase {
 public:
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

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

  ERICore(const Atomic& integral_calculator, const SturmintSystem& system,
          IntegralType type)
        : IntegralCoreBase(system, {IntegralCollection::id, type}),
          m_integral_calculator(integral_calculator) {
    assert_dbg(type == IntegralType::exchange || type == IntegralType::coulomb,
               krims::ExcInternalError());
  }

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap                                     */
  virtual void update(const krims::GenMap& map) override;

  /** \brief Clone the expression */
  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new ERICore(*this));
  }

 private:
  const Atomic& m_integral_calculator;
};

}  // namespace cs_reference
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
