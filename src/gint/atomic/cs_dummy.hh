#pragma once
#include <linalgwrap/Exceptions.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>
#include <sturmint/atomic/data/cs_dummy.hh>
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <linalgwrap/ParameterMap.hh>
#include "../Integral.hh"
#include "../IntegralCoreBase.hh"

namespace gint {
namespace atomic {
namespace cs_dummy {

  using namespace sturmint::atomic::cs_dummy;

  template <typename StoredMatrix> class OverlapIntegralCore;
  template <typename StoredMatrix> class NuclearAttractionIntegralCore;
  template <typename StoredMatrix> class KineticIntegralCore;
  template <typename StoredMatrix> class ERICore;
  
  //TODO: This has OrbitalType COMPLEX_ATOMIC. Decide how to multiplex wrt OrbitalType later.
  template <typename StoredMatrix> class IntegralCollection {
  public:
    typedef StoredMatrix                 stored_matrix_type;
    typedef Integral<stored_matrix_type> integral_matrix_type;
    typedef sturmint::real_t             real_type;
    
    static std::string id, name;
    
    real_type k_exponent, Z_charge;
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
    Z_charge{parameters.at<double>("Z_charge")},
    n_max{parameters.at<int>("n_max")},
    l_max{parameters.at<int>("l_max")},
    integral_calculator{n_max,l_max}
  {
    
  }

  template <typename StoredMatrix>
  Integral<StoredMatrix >IntegralCollection<StoredMatrix>::operator()(const std::string& integral_name)
  {
    if(integral_name == "nuclear_attraction") return Integral<StoredMatrix>{make_unique<NuclearAttractionIntegralCore<stored_matrix_type>>(integral_calculator,k_exponent,Z_charge)};
    if(integral_name == "overlap")            return Integral<StoredMatrix>{make_unique<OverlapIntegralCore<stored_matrix_type>>(integral_calculator)};     
    if(integral_name == "kinetic")            return Integral<StoredMatrix>{make_unique<KineticIntegralCore<stored_matrix_type>>(integral_calculator,k_exponent)};
    if(integral_name == "coulomb")            return Integral<StoredMatrix>{make_unique<ERICore<stored_matrix_type>>(integral_calculator,false,k_exponent)};  
    if(integral_name == "exchange")           return Integral<StoredMatrix>{make_unique<ERICore<stored_matrix_type>>(integral_calculator,true,k_exponent)};

    assert_dbg(false, linalgwrap::ExcNotImplemented());
    return Integral<StoredMatrix>(nullptr);
  }

  // ----------------------------------------------------------------------
  //			    INTEGRAL CORES
  // ----------------------------------------------------------------------
template <typename StoredMatrix>
class NuclearAttractionIntegralCore: public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix>    base_type;
  typedef StoredMatrix             stored_matrix_type;
  typedef typename base_type::size_type     size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t                  real_type;

  const real_type k, Z;
  
  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;
    
    assert_size(n_cols(),X.n_rows());

    stored_matrix_type AX(n_rows(), X.n_cols(),false);

    // TODO: n,l,m -> m,l,n
    for(size_t i=0;i<n_rows();i++){
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n;
      
      for(size_type j=0;j<AX.n_cols();j++) // TODO: make this work on vectors instead.
	AX(i,j) = (-Z*k/n)*X(i,j);
    }

    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return k*m_integral_calculator.nuclear_attraction(mui,muj);
  }  


  NuclearAttractionIntegralCore(const Atomic& integral_calculator, real_type k, real_type Z) : k(k), Z(Z), m_integral_calculator(integral_calculator) { }

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }
  
  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new NuclearAttractionIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id()   const override { return "atomic/cs_dummy/nuclear_attraction"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Nuclear attraction operator"; }

  
private:
  const Atomic&   m_integral_calculator;
};

  
template <typename StoredMatrix>
class OverlapIntegralCore: public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix>    base_type;
  typedef StoredMatrix             stored_matrix_type;
  typedef typename base_type::size_type     size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t real_type;

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

      for(size_type j=0;j<AX.n_cols();j++) // TODO: make this work on vectors instead.
	AX(i,j) = ((n>1)?     m_nnl*X(i_nminus,j) : 0)
	        +             X(i,j)
   	        + ((n<nmax)?  m_nnl*X(i_nplus,j) : 0);
    }

    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

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
  std::string name() const override { return "Overlap operator"; }

private:
  const Atomic& m_integral_calculator;
};


template <typename StoredMatrix>
class KineticIntegralCore: public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix>    base_type;
  typedef StoredMatrix             stored_matrix_type;
  typedef typename base_type::size_type     size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t                  real_type;
  
  real_type k;		// k-exponent
  
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
      double T_nnl = -0.5L*m_integral_calculator.overlap(n,n+1,l); // Kinetic is symmetric in n1,n2

      for(size_type j=0;j<AX.n_cols();j++) { 
	AX(i,j) = ( (n>1)?    T_nnl*X(i_nminus,j) : 0 )
	        +    0.5L*X(i,j)
	        + ( (n<nmax)? T_nnl*X(i_nplus,j) : 0 );

	AX(i,j) *= k*k;
      }

    }

    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t
      mui = nlmbasis::quantum_numbers_from_index(row),
      muj = nlmbasis::quantum_numbers_from_index(col);

    return k*k*m_integral_calculator.kinetic(mui,muj);
  }  

  KineticIntegralCore(const Atomic& integral_calculator, real_type k)
    : k(k), m_integral_calculator(integral_calculator) { }

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }
  
  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new KineticIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id()   const override { return "atomic/cs_dummy/kinetic"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Kinetic energy operator"; }

private:
  const Atomic& m_integral_calculator;
};

template <typename StoredMatrix>
class ERICore: public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix>    base_type;
  typedef StoredMatrix             stored_matrix_type;
  typedef typename base_type::size_type     size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t                  real_type;

  bool exchange;		 // Is this exchange or Coulomb operator?
  real_type k;		         // Exponent scale
  std::shared_ptr<stored_matrix_type>  coefficients_occupied;
  
  
  /** \brief Multiplication with a stored matrix */
  // J_{b1,q} = J_{b1,b2} X_{b2,q} = J_{b1,b2,b3,b4} X_{b2,q} Cocc_{b3,p} Cocc_{b4,p} = J_{b1,b2,b3,b4} X_{b2,q} D_{b3,b4}
  // K_{b1,q} = K_{b1,b2} X_{b2,q} = J_{b1,b2,b3,b4} X_{b3,q} Cocc_{b2,p} Cocc_{b4,p} = J_{b1,b2,b3,b4} X_{b2,q} D_{b3,b4}
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    assert_size(n_cols(),X.n_rows());
      
    stored_matrix_type AX(n_rows(), X.n_cols(),true);
    const int l_max = m_integral_calculator.lmax1;
    const int n_max = m_integral_calculator.nmax1;
    const stored_matrix_type &Cocc(*coefficients_occupied);

    cout << "\nCalculating "<<(exchange?"exchange":"coulomb")<<" integral application.\n";
    for(size_t b1=0;b1<n_rows();b1++){ // TODO: Mål, om det bliver hurtigere / nemmere at OpenMP'e, af at iterere over b1q samlet.
      nlm_t nlm1 = nlmbasis::quantum_numbers_from_index(b1);
      int8_t /*n1 = nlm1.n, */l1 = nlm1.l, m1 = nlm1.m;
	
      for(size_t q=0;q<X.n_cols();q++){
	real_type sum = 0;
	  
	for(size_t b2=0;b2<X.n_rows();b2++){
	  nlm_t nlm2 = nlmbasis::quantum_numbers_from_index(b2);
	  int8_t /*n2 = nlm2.n, */l2 = nlm2.l, m2 = nlm2.m;

	  for(size_t b3=0;b3<n_rows();b3++){
	    nlm_t nlm3 = nlmbasis::quantum_numbers_from_index(b3);
	    int8_t n3 = nlm3.n, l3 = nlm3.l, m3 = nlm3.m;
	    cout << "b3 = " << b3 << " = index(" << vector<int>{n3,l3,m3} << ")\n";

	    // TODO: Exchange? Da!	    
	    int8_t m4 = m3-m2+m1;
	    int l_parity = (l1+l2)&1, m_parity = (m1+m2)&1;
	    int l_min    = max(l_parity, ::abs(m4)+((m_parity+l_parity)&1)); // TODO: Check!!
	    cout << "\t{l_parity,m_parity,l_min} = " << vector<int>{l_parity,m_parity,l_min} << endl;
	    
	    int B2 = exchange? b3 : b2; // Swap b2 and b3 if computing exchange 
	    int B3 = exchange? b2 : b3;
	    real_type X_bq = X(B2,q);
	    
	    for(int8_t l4=l_min;l4<=l_max;l4+=2){
	      //		int8_t lmin = max(abs(m2-m1)+(m_parity-l_parity), max(abs(l1-l2),abs(l3-l4))); // l >= |m2-m1| && l >= |l1-l2| && l >= |l3-l4| && (-1)^l = (-1)^(l1+l2)
	      //		int8_t lmax = min(l1+l2,l3+l4);

	      for(int8_t n4=l4+1;n4<=n_max;n4++){ // TODO: Fixed number of n's (same-length dot product for all (l,m)-combinations). Same (l,m), different n's should be contiguous in memory.
		// TODO: Ide: Shell-pairs grupperet efter (m2-m1) og l-paritet. Mål: contiguous dot-product.
		// TODO: Alle atomic/cs implementationer er fuldstændig ens bortset fra repulsion-operationerne.
		int b4 = nlmbasis::index(nlm_t{n4,l4,m4});
		cout << "\tb4 = " << b4 << " = index(" << vector<int>{n4,l4,m4} << ")\n";
		  
		for(size_t p=0;p<Cocc.n_cols();p++)
		  sum += sturmint::atomic::cs_dummy::repulsion[b4+norb*(b3+norb*(b2+norb*b1))]*Cocc(B3,p)*Cocc(b4,p)*X_bq;
	      }
	    }
	  }
	}

	AX(b1,q) = sum;
      }
    }
	
    return AX;
  }

  /** \brief return an element of the matrix    */
  // J_{b1,b2} = J_{b1,b2,b3,b4} Cocc_{b3,p} Cocc_{b4,p} = J_{b1,b2,b3,b4} D_{b3,b4}
  scalar_type operator()(size_type b1, size_type b2) const override {
    using namespace sturmint::atomic;
    using namespace sturmint::orbital_index;
    
    const int l_max = m_integral_calculator.lmax1;
    const int n_max = m_integral_calculator.nmax1;
    const stored_matrix_type &Cocc(*coefficients_occupied);
      
    nlm_t nlm1 = nlmbasis::quantum_numbers_from_index(b1);
    nlm_t nlm2 = nlmbasis::quantum_numbers_from_index(b2);
    int8_t n1 = nlm1.n, l1 = nlm1.l, m1 = nlm1.m;
    int8_t n2 = nlm2.n, l2 = nlm2.l, m2 = nlm2.m;

    real_type sum = 0;

    cout << "\nCalculating "<<(exchange?"exchange":"coulomb")<<" integral element ("<<b1<<","<<b2<<").\n"
	 << "{n1,l1,m1} = " << vector<int>{n1,l1,m1} << "; {n2,l2,m2} = " << vector<int>{n2,l2,m2} << ";\n";
    
    for(size_t b3=0;b3<n_rows();b3++){
      nlm_t nlm3 = nlmbasis::quantum_numbers_from_index(b3);
      int8_t n3 = nlm3.n, l3 = nlm3.l, m3 = nlm3.m;

      int8_t m4 = m3-m2+m1;
      int l_parity = (l1+l2)&1, m_parity = (m1+m2)&1;
      int l_min    = max(l_parity, ::abs(m4)+((m_parity+l_parity)&1)); // TODO: Check!!
      cout << "b3 = " << b3 << " = index(" << vector<int>{n3,l3,m3} << ")\n";      
      cout << "\t{l_parity,m_parity,l_min} = " << vector<int>{l_parity,m_parity,l_min} << endl;

      //      int B2 = exchange? b3 : b2; // Swap b2 and b3 if computing exchange // TODO: Check
      int B3 = exchange? b2 : b3;

      for(int8_t l4=l_min;l4<=l_max;l4+=2){
	//		int8_t lmin = max(abs(m2-m1)+(m_parity-l_parity), max(abs(l1-l2),abs(l3-l4))); // l >= |m2-m1| && l >= |l1-l2| && l >= |l3-l4| && (-1)^l = (-1)^(l1+l2)
	//		int8_t lmax = min(l1+l2,l3+l4);

	// TODO: Fixed number of n's (same-length dot product for all
	// (l,m)-combinations). Same (l,m), different n's should be
	// contiguous in memory.
	for(int8_t n4=l4+1;n4<=n_max;n4++){ 
	  int b4 = nlmbasis::index(nlm_t{n4,l4,m4});
	  cout << "\tb4 = " << b4 << " = index(" << vector<int>{n4,l4,m4} << ")\n";

	  for(size_t p=0;p<Cocc.n_cols();p++)
	    sum += cs_dummy::repulsion[b4+norb*(b3+norb*(b2+norb*b1))]*Cocc(B3,p)*Cocc(b4,p);
	}
      }
    }
    return sum;
  }
	
  ERICore(const Atomic& integral_calculator, bool exchange, real_type k) : exchange(exchange), k(k), m_integral_calculator(integral_calculator) { }

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap                                     */
  virtual void update(const linalgwrap::ParameterMap& p) {
    // TODO: Don't copy coefficients, copy a shared pointer instead
    // coefficients_occupied = p.at<stored_matrix_type>("coefficients_occupied");
    coefficients_occupied = make_shared<stored_matrix_type>(p.at<stored_matrix_type>("coefficients_occupied"));
  }

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }
  
  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new ERICore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id()   const override { return std::string("atomic/cs_dummy/ERI_")+(exchange?"K":"J"); }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return std::string("Electron Repulsion Integrals, ")+(exchange?"Exchange":"Coulomb")+" operator"; }

private:
  const Atomic& m_integral_calculator;
};
  
  
} // namespace cs_dummy
} // namespace atomic
} // namespace gint
