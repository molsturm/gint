#pragma once
#include "gint/IntegralCollectionBase.hh"
#include "nlm_order.hh"
#include <sturmint/atomic/cs_naive/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_naive {
using namespace nlm_order;
using sturmint::atomic::cs_naive::Atomic;

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
   */
  IntegralCollection(const krims::GenMap& parameters);

  using base_type::lookup_integral;
  integral_matrix_type lookup_integral(IntegralType type) const override;

  /** Obtain the id string of the collection / basis type */
  const std::string& basis_id() const override { return id; }

  /** Obtain the friendly name of the collection / basis type */
  std::string basis_name() const override {
    return "Naive implementation of atomic Coulomb Sturmians";
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

class ERICore final : public nlm_order::ERICore<Atomic> {
 public:
  using nlm_order::ERICore<Atomic>::ERICore;

  /** \brief Multiplication with a stored matrix */
  // J_{aq} = J_{ab} X_{bq} = J_{abcd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{abcd} X_{bq} D_{cd}
  // K_{aq} = K_{ab} X_{bq} = J_{acbd} X_{bq} Cocc_{cp} Cocc_{dp} = J_{acbd} X_{bq} D_{cd}
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** \brief Clone the expression */
  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new ERICore(*this));
  }
};

}  // namespace cs_naive
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
