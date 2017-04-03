#pragma once
#ifdef GINT_HAVE_LIBINT

#include <functional>
#include <krims/SubscriptionPointer.hh>
#include <libint2.hpp>

#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/Structure.hh"
#include "gint/config.hh"

namespace gint {
namespace gaussian {
namespace libint {

// In this namespace all things are real:
using namespace real_valued;

class OverlapIntegralCore;
class NuclearAttractionIntegralCore;
class KineticIntegralCore;
class ERICore;

/** Initialises the global libint state and manages it */
class LibintGlobalInit : public krims::Subscribable {
 public:
  LibintGlobalInit() { libint2::initialize(); }
  ~LibintGlobalInit() { libint2::finalize(); }
  LibintGlobalInit(LibintGlobalInit&&) = default;
  LibintGlobalInit& operator=(LibintGlobalInit&&) = default;
};

/** Object representing a molecular structure and basis in a way usable to libint.
 *
 * In a way this is the system (physical and model) we want to compute integrals for.
 **/
class LibintSystem : public krims::Subscribable {
 public:
  LibintSystem() : m_structure_ptr("LibintSystem") {}
  LibintSystem(const std::string& basisset_name,
               krims::SubscriptionPointer<const Structure> structure_ptr);

  /** Access to the libint basis structure */
  const libint2::BasisSet& basis() const { return m_basis; }

  /** Access the molecular structure, which was used to build this basis */
  const Structure& structure() const { return *m_structure_ptr; }

  /** Access the list of point charges used by libint */
  const std::vector<std::pair<double, std::array<double, 3>>>& point_charges() const {
    return m_point_charges;
  }

  /** Access to the number of shells (1s, 2s, 2p, ...)*/
  size_t n_shells() const { return m_basis.size(); }

  /** Access to the total number of basis functions */
  size_t n_bas() const { return static_cast<size_t>(m_basis.nbf()); }

 private:
  libint2::BasisSet m_basis;
  krims::SubscriptionPointer<const Structure> m_structure_ptr;
  std::vector<std::pair<double, std::array<double, 3>>> m_point_charges;
};

/** Information about a particular shell */
struct LibintShell {
  //! Shell index
  size_t index;

  //! Number of basis functions in this shell
  size_t n_bfct;

  //! Index of the first basis function of this shell
  size_t first_bfct;
};

/** IntegralCollection for the libint integral objects */
class IntegralCollection final : public IntegralCollectionBase<stored_matrix_type> {
 public:
  typedef IntegralCollectionBase<stored_matrix_type> base_type;

  const static std::string id;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - basis_set  (std::string)  The string describing the
   *                                    Gaussian basis set to use.
   *   - structure (gint::Structure)  The structure of the molecule.
   */
  explicit IntegralCollection(const krims::GenMap& parameters);

  using base_type::lookup_integral;
  integral_matrix_type lookup_integral(IntegralType type) const override;

  /** Obtain the id string of the collection / basis type */
  const std::string& basis_id() const override { return id; }

  /** Obtain the friendly name of the collection / basis type */
  std::string basis_name() const override { return "Gaussian integrals from libint2"; }

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::unique_ptr<base_type> create(const krims::GenMap& parameters) {
    return krims::make_unique<IntegralCollection>(parameters);
  }

 private:
  LibintSystem m_system;
  LibintGlobalInit m_global;
};

//
// Integral Cores
//

/** Base class for libint integral cores */
class LibintIntegralCoreBase : public IntegralCoreBase<stored_matrix_type> {
 public:
  typedef IntegralCoreBase<stored_matrix_type> base_type;
  typedef typename base_type::scalar_type scalar_type;

  /** Constructor
   *
   * \param system  The molecular system and basis this core deals with
   * \param global  The global libint resources to acquire
   */
  LibintIntegralCoreBase(const LibintSystem& system, const LibintGlobalInit& global)
        : m_system_ptr("LibIntOneElectronIntegralCore", system),
          m_global_ptr("LibIntOneElectronIntegralCore", global) {}

  size_t n_rows() const override { return m_system_ptr->n_bas(); }
  size_t n_cols() const override { return m_system_ptr->n_bas(); }

  scalar_type operator()(size_t row, size_t col) const override;

  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A = 1, const scalar_type c_M = 0) const override;

  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

 protected:
  /** Return the system object */
  const LibintSystem& system() const { return *m_system_ptr; }

  /** Type of a computational kernel for the compute function below.
   *
   * For each pair of shell indices (e.g. 1s, 2s, 2p, ...) the kernel function
   * is called once.
   */
  typedef std::function<void(LibintShell, LibintShell, const scalar_type*)> kernel_type;

  // TODO Calling a lambda is faster than calling via a std::function object.
  //      Try to think of a way to achieve this.

  //@{
  /** Call a kernel function for each pair of shells of the basis set.
   *
   * Optionally one can supply ranges of *basis* indices in rows and columns.
   * The function will then only loop over shell indices which cover the requested
   * range of rows and columns.
   *
   * Note that the kernel may well be called with a larger amount of data present
   * than what was requested.
   *
   * \param row_shells   The range of row shell indices to work on (optional)
   * \param col_shells   The range of column shell indices to work on (optional)
   * \param kernel  The kernel function to excecute for each pair of shell indices
   *                This function should take 3 arguments:
   *                        - size_t            Row indices to compute data for
   *                        - size_t            Col indices to compute data for
   *                        - const scalar_type*  Buffer of the integral values.
   **/
  virtual void compute(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
                       kernel_type&& kernel) const = 0;

  void compute(kernel_type&& kernel) const {
    assert_dbg(n_rows() == n_cols(), krims::ExcInternalError());
    const auto full_range = krims::range(n_rows());
    compute(full_range, full_range, std::forward<kernel_type>(kernel));
  }
  //@}

 private:
  // Pointer to the molecular system info we use:
  krims::SubscriptionPointer<const LibintSystem> m_system_ptr;

  // Pointer to the global libint2 state we need
  krims::SubscriptionPointer<const LibintGlobalInit> m_global_ptr;
};

/** Class for one electron integral cores with libint */
class OneElecIntegralCore final : public LibintIntegralCoreBase {
 public:
  typedef LibintIntegralCoreBase base_type;
  typedef IntegralCoreBase<stored_matrix_type> core_base_type;

  /** Constructor
   *
   * \param system  The molecular system and basis this core deals with
   * \param op      The operator libint should compute in this core
   * \param global  The global libint resources to acquire
   */
  OneElecIntegralCore(libint2::Operator op, const LibintSystem& system,
                      const LibintGlobalInit& global)
        : base_type(system, global), m_operator(op) {}

  std::unique_ptr<core_base_type> clone() const override {
    return std::unique_ptr<core_base_type>(new OneElecIntegralCore(*this));
  }

  IntegralIdentifier id() const override {
    switch (m_operator) {
      case libint2::Operator::overlap:
        return {IntegralCollection::id, IntegralType::overlap};
      case libint2::Operator::kinetic:
        return {IntegralCollection::id, IntegralType::kinetic};
      case libint2::Operator::nuclear:
        return {IntegralCollection::id, IntegralType::nuclear_attraction};
      default:
        // This is not a one electron integral operator known to us atm.
        assert_throw(false, krims::ExcInternalError());
        return {IntegralCollection::id, IntegralType::overlap};
    }
  }

 protected:
  virtual void compute(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
                       typename base_type::kernel_type&& kernel) const override;

 private:
  // Integral operator this core computes.
  libint2::Operator m_operator;
};

class ERICore final : public LibintIntegralCoreBase {
 public:
  typedef LibintIntegralCoreBase base_type;
  typedef IntegralCoreBase<stored_matrix_type> core_base_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  ERICore(IntegralType type, const LibintSystem& system, const LibintGlobalInit& global)
        : base_type(system, global), m_type(type) {
    assert_dbg(type == IntegralType::exchange || type == IntegralType::coulomb,
               krims::ExcInternalError());
  }

  std::unique_ptr<core_base_type> clone() const override {
    return std::unique_ptr<core_base_type>(new ERICore(*this));
  }

  IntegralIdentifier id() const override { return {IntegralCollection::id, m_type}; }

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap
   **/
  void update(const krims::GenMap& map) override;

 protected:
  void compute(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
               kernel_type&& kernel) const override;

  IntegralType m_type;
};

}  // namespace libint
}  // namespace gaussian
}  // namespace gint
#endif  // GINT_HAVE_LIBINT
