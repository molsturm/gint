#pragma once
#include "gint/IntegralCollectionBase.hh"
#include "nlm_order.hh"
#include <sturmint/atomic/cs_reference_pc/cs_atomic.hh>

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_reference_pc {
using namespace nlm_order;
using sturmint::atomic::cs_reference_pc::Atomic;

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

class ERICore : public ERICoreBase {
 public:
  /** \brief return an element of the matrix    */
  // J_{ab} = J_{abcd} Cocc_{cp} Cocc_{dp} = J_{abcd} P_{cd}
  // K_{ab} = J_{cbad} Cocc_{cp} Cocc_{dp} = J_{acbd} P_{cd}
  scalar_type operator()(size_t a, size_t b) const override;

  ERICore(const Atomic& integral_calculator, const SturmintSystem& system,
          IntegralType type)
        : ERICoreBase(system, {IntegralCollection::id, type}),
          m_integral_calculator(integral_calculator) {
    assert_dbg(type == IntegralType::exchange || type == IntegralType::coulomb,
               krims::ExcInternalError());
  }

  /** \brief Clone the expression */
  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new ERICore(*this));
  }

 private:
  const Atomic& m_integral_calculator;
};

}  // namespace cs_reference_pc
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
