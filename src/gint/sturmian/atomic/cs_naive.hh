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

}  // namespace cs_naive
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
