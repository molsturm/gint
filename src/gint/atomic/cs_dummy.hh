#pragma once
#include "../Integral.hh"
#include "../IntegralCoreBase.hh"
#include <linalgwrap/Exceptions.hh>
#include <linalgwrap/ParameterMap.hh>
#include <sturmint/atomic/cs_dummy/cs_atomic.hh>
#include <sturmint/harmonic/OrbitalIndex.hh>

namespace gint {
namespace atomic {
namespace cs_dummy {
using namespace sturmint::atomic::cs_dummy;

// Forward declare integral cores
template <typename StoredMatrix>
class OverlapIntegralCore;
template <typename StoredMatrix>
class NuclearAttractionIntegralCore;
template <typename StoredMatrix>
class KineticIntegralCore;
template <typename StoredMatrix>
class JIntegralCore;
template <typename StoredMatrix>
class KIntegralCore;

// TODO: This has OrbitalType COMPLEX_ATOMIC. Decide how to multiplex wrt
// OrbitalType later.
template <typename StoredMatrix>
class IntegralCollection {
public:
  typedef StoredMatrix stored_matrix_type;
  typedef Integral<stored_matrix_type> integral_matrix_type;
  // TODO do this typedef differently (e.g. by having a common real type)
  typedef sturmint::real_t real_type;

  //! Id string of this type of integrals
  static std::string id;

  //! Friendly name of this type of integrals
  static std::string name;

  //! Exponent of all Coulomb sturmians
  real_type k_exponent;

  //! Maximum principle quantum number used by the basis
  int n_max;

  //! Maximum angular momentum quantum number used by the basis
  int l_max;

  //! (Fractional) nuclear charge of the system
  real_type Z_charge;

  //! The object building the matrix elements.
  Atomic integral_calculator;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - Z_charge (double): The nuclear change of the system
   *   - n_max (int): The maximal principle quantum number
   *   - l_max (int): Maximal angular momentum qn
   */
  IntegralCollection(const linalgwrap::ParameterMap& parameters);

  /** Lookup an integral by its identifier string */
  integral_matrix_type operator()(const std::string& integral_id);
};

// ----------------------------------------------------------------------
//			    IMPLEMENTATION
// ----------------------------------------------------------------------
template <typename StoredMatrix>
std::string IntegralCollection<StoredMatrix>::id = "cs_atomic/dummy";
template <typename StoredMatrix>
std::string IntegralCollection<StoredMatrix>::name =
      "Dummy implementation of atomic Coulomb Sturmians";

template <typename StoredMatrix>
IntegralCollection<StoredMatrix>::IntegralCollection(
      const linalgwrap::ParameterMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        n_max{parameters.at<int>("n_max")},
        l_max{parameters.at<int>("l_max")},
        Z_charge{parameters.at<double>("Z_charge")},
        integral_calculator{n_max, l_max} {
  assert_throw(n_max > 0,
               ExcInvalidIntegralParameters(
                     "Maximum principle quantum number (" +
                     std::to_string(n_max) + ") needs to be greater 0."));
  assert_throw(l_max >= 0, ExcInvalidIntegralParameters(
                                 "Maximum angular momentum quantum number (" +
                                 std::to_string(l_max) +
                                 ") needs to be greater or equal 0."));
  assert_throw(l_max < n_max,
               ExcInvalidIntegralParameters(
                     "Maximum angular momentum (" + std::to_string(l_max) +
                     ") needs to be smaller than maximum principle "
                     "quantum number (" +
                     std::to_string(n_max) + ")."));
}

template <typename StoredMatrix>
Integral<StoredMatrix> IntegralCollection<StoredMatrix>::operator()(
      const std::string& integral_name) {
  if (integral_name == "nuclear_attraction")
    return make_integral<NuclearAttractionIntegralCore<stored_matrix_type>>(
          integral_calculator, k_exponent, Z_charge);
  if (integral_name == "overlap")
    return make_integral<OverlapIntegralCore<stored_matrix_type>>(
          integral_calculator);
  if (integral_name == "kinetic")
    return make_integral<KineticIntegralCore<stored_matrix_type>>(
          integral_calculator, k_exponent);
  if (integral_name == "coulomb")
    return make_integral<JIntegralCore<stored_matrix_type>>(integral_calculator,
                                                            k_exponent);
  if (integral_name == "exchange")
    return make_integral<KIntegralCore<stored_matrix_type>>(integral_calculator,
                                                            k_exponent);

  assert_dbg(false, linalgwrap::ExcNotImplemented());
  return Integral<StoredMatrix>(nullptr);
}

// ----------------------------------------------------------------------
//                           INTEGRAL CORES
// ----------------------------------------------------------------------

template <typename StoredMatrix>
class NuclearAttractionIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t real_type;

  const real_type k, Z;

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;

    assert_size(n_cols(), X.n_rows());

    stored_matrix_type AX(n_rows(), X.n_cols(), false);

    // TODO: n,l,m -> m,l,n
    for (size_t i = 0; i < n_rows(); i++) {
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n;

      for (size_type j = 0; j < AX.n_cols();
           j++)  // TODO: make this work on vectors instead.
                 // Also TODO: Make the memory layout more friendly to operator
                 // application on many vectors.
        AX(i, j) = (-Z * k / n) * X(i, j);
    }

    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return k * m_integral_calculator.nuclear_attraction(mui, muj);
  }

  NuclearAttractionIntegralCore(const Atomic& integral_calculator, real_type k,
                                real_type Z)
        : k(k), Z(Z), m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }

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

template <typename StoredMatrix>
class OverlapIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t real_type;

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;

    assert_size(n_cols(), X.n_rows());

    stored_matrix_type AX(n_rows(), X.n_cols(), false);
    const int nmax = m_integral_calculator.nmax1;

    // TODO: n,l,m -> m,l,n
    for (size_t i = 0; i < n_rows(); i++) {
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n, l = nlm.l, m = nlm.m;

      size_t i_nminus = nlmbasis::index(nlm_t(n - 1, l, m)),
             i_nplus = nlmbasis::index(nlm_t(n + 1, l, m));

      double m_nnl = m_integral_calculator.overlap(
            n, n + 1, l);  // Overlap is symmetric in n,np

      for (size_type j = 0; j < AX.n_cols();
           j++)  // TODO: make this work on vectors instead.
                 // Also TODO: Make the memory layout more friendly to operator
                 // application on many vectors.
        AX(i, j) = ((n > 1) ? m_nnl * X(i_nminus, j) : 0) + X(i, j) +
                   ((n < nmax) ? m_nnl * X(i_nplus, j) : 0);
    }

    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return m_integral_calculator.overlap(mui, muj);
  }

  OverlapIntegralCore(const Atomic& integral_calculator)
        : m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }

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

template <typename StoredMatrix>
class KineticIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t real_type;

  real_type k;  // k-exponent

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    typedef nlmbasis::quantum_numbers_t nlm_t;

    assert_size(n_cols(), X.n_rows());

    stored_matrix_type AX(n_rows(), X.n_cols(), false);
    const int nmax = m_integral_calculator.nmax1;

    // TODO: n,l,m -> m,l,n
    for (size_t i = 0; i < n_rows(); i++) {
      nlm_t nlm = nlmbasis::quantum_numbers_from_index(i);
      int8_t n = nlm.n, l = nlm.l, m = nlm.m;

      size_t i_nminus = nlmbasis::index(nlm_t(n - 1, l, m)),
             i_nplus = nlmbasis::index(nlm_t(n + 1, l, m));

      double T_nnl = -0.5L *
                     m_integral_calculator.overlap(
                           n, n + 1, l);  // Kinetic is symmetric in n,np

      for (size_type j = 0; j < AX.n_cols();
           j++) {  // TODO: make this work on vectors instead.
        AX(i, j) = ((n > 1) ? T_nnl * X(i_nminus, j) : 0) + 0.5L * X(i, j) +
                   ((n < nmax) ? T_nnl * X(i_nplus, j) : 0);

        AX(i, j) *= k * k;
      }
    }

    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    using sturmint::orbital_index::nlmbasis;
    const nlm_t mui = nlmbasis::quantum_numbers_from_index(row),
                muj = nlmbasis::quantum_numbers_from_index(col);

    return k * k * m_integral_calculator.kinetic(mui, muj);
  }

  KineticIntegralCore(const Atomic& integral_calculator, real_type k)
        : k(k), m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }

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

template <typename StoredMatrix>
class JIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t real_type;

  real_type k;  // k-exponent

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    assert_size(n_cols(), X.n_rows());
    stored_matrix_type AX(n_rows(), X.n_cols(), true);
    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    return 0;
  }

  JIntegralCore(const Atomic& integral_calculator, real_type k)
        : k(k), m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new JIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id() const override { return "atomic/cs_dummy/J"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Coulomb operator J"; }

private:
  const Atomic& m_integral_calculator;
};

template <typename StoredMatrix>
class KIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef sturmint::real_t real_type;

  real_type k;  // k-exponent

  /** \brief Multiplication with a stored matrix */
  stored_matrix_type operator*(const stored_matrix_type& X) const override {
    using namespace sturmint::orbital_index;
    assert_size(n_cols(), X.n_rows());
    stored_matrix_type AX(n_rows(), X.n_cols(), true);
    return AX;
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override {
    return 0;
  }

  KIntegralCore(const Atomic& integral_calculator, real_type k)
        : k(k), m_integral_calculator(integral_calculator) {}

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_integral_calculator.n_bas(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_integral_calculator.n_bas(); }

  /** \brief Clone the expression */
  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new KIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  std::string id() const override { return "atomic/cs_dummy/K"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Exchange operator K"; }

private:
  const Atomic& m_integral_calculator;
};

}  // namespace cs_dummy
}  // namespace atomic
}  // namespace gint
