#pragma once

namespace gint {
namespace atomic {
namespace cs_dummy {
  using namespace sturmint::atomic::cs_dummy;
  
  //TODO: This has OrbitalType COMPLEX_ATOMIC. Decide how to multiplex wrt OrbitalType later.
  templace <typename StoredMatrix> class DummyCollection {
  public:
    typedef StoredMatrix                 stored_matrix_type;
    typedef Integral<stored_matrix_type> integral_matrix_type;
    
    static std::string id, name;
    
    IntegralMatrixType nuclear_attraction_matrix, overlap_matrix,  kinetic_energy_matrix, J_matrix, K_matrix;

    double k_exponent;
    int    n_max, l_max;
    Atomic integral_calculator;

    DummyIntegrals(const linalgwrap::ParameterMap& parameters);
    integral_matrix_type operator()(const std::string& integral_name);
  };

  // ----------------------------------------------------------------------
  //			    IMPLEMENTATION
  // ----------------------------------------------------------------------
  template <typename StoredMatrix>
  std::string DummyIntegrals::id   = "cs_atomic/dummy";
  std::string DummyIntegrals::name = "Dummy implementation of atomic Coulomb Sturmians";
  
  template <typename StoredMatrix> 
  DummyIntegrals::DummyIntegrals(const linalgwrap::ParameterMap& parameters) :
    k_exponent{parameters.get<double>("k_exponent")},
    n_max{parameters.get<int>("n_max")},
    l_max{parameters.get<int>("l_max")},
    integral_calculator{n_max,l_max},
    overlap_matrix{make_unique<OverlapIntegralCore<stored_matrix_type>>(integral_calculator)}
  {
    
  }

  template <typename StoredMatrix>
  Integral<StoredMatrix> DummyIntegrals::operator()(const std::string& integral_name)
  {
    if(integral_name == "nuclear_attraction") return nuclear_attraction_matrix;
    if(integral_name == "overlap")            return overlap_matrix;
    if(integral_name == "kinetic")            return kinetic_energy_matrix;
    if(integral_name == "coulomb")            return J_matrix;
    if(integral_name == "exchange")           return K_matrix;

    assert_dbg(false, ExcNotImplemented());
  }

} // namespace cs_dummy
} // namespace atomic
} // namespace gint
