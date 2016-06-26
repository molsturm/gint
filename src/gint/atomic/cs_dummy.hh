#pragma once
#include <sturmint/harmonic/OrbitalIndex.hh>
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <linalgwrap/ParameterMap.hh>
#include "../Integral.hh"
#include "../IntegralCoreBase.hh"


namespace gint {
namespace atomic {
namespace cs_dummy {

  using namespace sturmint::atomic::cs_dummy;

  template <typename StoredMatrix> class OverlapIntegralCore;

  
  //TODO: This has OrbitalType COMPLEX_ATOMIC. Decide how to multiplex wrt OrbitalType later.
  template <typename StoredMatrix> class IntegralCollection {
  public:
    typedef StoredMatrix                 stored_matrix_type;
    typedef Integral<stored_matrix_type> integral_matrix_type;
    
    static std::string id, name;
    
    integral_matrix_type nuclear_attraction_matrix, overlap_matrix,  kinetic_energy_matrix, J_matrix, K_matrix;

    double k_exponent;
    int    n_max, l_max;
    Atomic integral_calculator;

    IntegralCollection(const linalgwrap::ParameterMap& parameters);
    integral_matrix_type operator()(const std::string& integral_name);
  };

  // ----------------------------------------------------------------------
  //			    IMPLEMENTATION
  // ----------------------------------------------------------------------
  template <typename StoredMatrix>
  std::string IntegralCollection<StoredMatrix>::id   = "cs_atomic/dummy";
  template <typename StoredMatrix>
  std::string IntegralCollection<StoredMatrix>::name = "Dummy implementation of atomic Coulomb Sturmians";
  
  template <typename StoredMatrix> 
  IntegralCollection<StoredMatrix>::IntegralCollection(const linalgwrap::ParameterMap& parameters) :
    k_exponent{parameters.at<double>("k_exponent")},
    n_max{parameters.at<int>("n_max")},
    l_max{parameters.at<int>("l_max")},
    integral_calculator{n_max,l_max},
    overlap_matrix{make_unique<OverlapIntegralCore<stored_matrix_type>>(integral_calculator)}
  {
    
  }

  template <typename StoredMatrix>
  Integral<StoredMatrix> IntegralCollection<StoredMatrix>::operator()(const std::string& integral_name)
  {
    if(integral_name == "nuclear_attraction") return nuclear_attraction_matrix;
    if(integral_name == "overlap")            return overlap_matrix;
    if(integral_name == "kinetic")            return kinetic_energy_matrix;
    if(integral_name == "coulomb")            return J_matrix;
    if(integral_name == "exchange")           return K_matrix;

    assert_dbg(false, ExcNotImplemented());
  }

  // ----------------------------------------------------------------------
  //			    INTEGRAL CORES
  // ----------------------------------------------------------------------

template <typename StoredMatrix>
class OverlapIntegralCore: public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix>    base_type;
  typedef StoredMatrix             stored_matrix_type;
  typedef typename base_type::size_type     size_type;
  typedef typename base_type::scalar_type scalar_type;

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;
    
    assert_size(n_cols(),X.n_rows());

    stored_matrix_type AX(n_rows(), X.n_cols(),false);
    const int nmax = m_integral_calculator.nmax1;

    // TODO: n,l,m -> m,l,n
    for(size_t i=0;i<n_rows();i++){
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n, l = nlm.l, m = nlm.m;
      
      size_t i_nminus = nlmbasis::index(nlm_t(n-1,l,m)), i_nplus = nlmbasis::index(nlm_t(n+1,l,m));

      double m_nnl = m_integral_calculator.overlap(n,n+1,l); // Overlap is symmetric in n,np

      for(int j=0;j<n_cols();j++) // TODO: make this work on vectors instead.
	                          // Also TODO: Make the memory layout more friendly to operator application on many vectors.
	AX(i,j) = (n>1)*m_nnl*X(i_nminus,j)  + X(i,j) + (n<nmax)*m_nnl*X(i_nplus,j);
    }

    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = quantum_numbers_from_index(row), muj = quantum_numbers_from_index(col);

    return m_integral_calculator.overlap(mui,muj);
  }  


  OverlapIntegralCore(const Atomic& integral_calculator) : m_integral_calculator(integral_calculator) { }

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }
  
  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new OverlapIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id()   const override { return "atomic/cs_dummy/overlap"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Overlap integral"; }

private:
  const Atomic& m_integral_calculator;
};

  
} // namespace cs_dummy
} // namespace atomic
} // namespace gint
