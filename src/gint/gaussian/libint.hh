#pragma once
#ifdef GINT_HAVE_LIBINT

#include <krims/SubscriptionPointer.hh>
#include <libint2.hpp>

#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/chemistry/Molecule.hh"
#include "gint/config.hh"

namespace gint {
namespace gaussian {
namespace libint {

// In this namespace all things are real:
typedef real_type scalar_type;
typedef real_stored_mtx_type stored_mtx_type;
typedef real_multivector_type multivector_type;
typedef const_real_multivector_type const_multivector_type;

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
               krims::SubscriptionPointer<const Molecule> structure_ptr);

  /** Access to the libint basis structure */
  const libint2::BasisSet& basis() const { return m_basis; }

  /** Access the molecular structure, which was used to build this basis */
  const Molecule& structure() const { return *m_structure_ptr; }

  /** Access the list of point charges used by libint */
  const std::vector<std::pair<double, std::array<double, 3>>>& point_charges() const {
    return m_point_charges;
  }

 private:
  libint2::BasisSet m_basis;
  krims::SubscriptionPointer<const Molecule> m_structure_ptr;
  std::vector<std::pair<double, std::array<double, 3>>> m_point_charges;
};

/** When initialised with a LibintSystem it computes and makes available some data which
 * is needed to work with libint shell pairs */
struct LibintShellData {
  /** Return the index of the first basis function of a particular shell */
  size_t first_bfct(size_t shell) const { return shell2bf[shell]; }

  /** Return the number of basis functions of a particular shell */
  size_t n_bfct(size_t shell) const { return m_system_ptr->basis()[shell].size(); }

  /** Return the shell range which covers the requested index range
   * in the basis functions.
   *
   * In other words if only some basis function indices are required,
   * this function returns the range of shells which need to be computed
   * in order to compute the values for all basis functions.
   * Note that it is very likely that extra values at the beginning and
   * the end will be computed as well.
   */
  krims::Range<size_t> shell_range(krims::Range<size_t> bfct_indices) const;

  /** The map shell index to index of first basis function of the shell
   *  Obtain the map shell index to basis function index
   *  E.g. shell2bf[0] = index of the first basis function in shell 0
   *       shell2bf[1] = index of the first basis function in shell 1
   */
  const std::vector<size_t> shell2bf;

  explicit LibintShellData(const LibintSystem& system)
        : shell2bf(system.basis().shell2bf()), m_system_ptr("LibintBasisData", system) {}

 private:
  // Pointer to the molecular system info we use:
  krims::SubscriptionPointer<const LibintSystem> m_system_ptr;
};

/** IntegralCollection for the libint integral objects */
class IntegralCollection final
      : public IntegralCollectionBase<OrbitalType::REAL_MOLECULAR> {
 public:
  typedef IntegralCollectionBase<OrbitalType::REAL_MOLECULAR> base_type;

  const static std::string id;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - basis_set  (std::string)  The string describing the
   *                                    Gaussian basis set to use.
   *   - structure (gint::Molecule)  The structure of the molecule.
   */
  explicit IntegralCollection(const krims::GenMap& parameters);

  /** Lookup an integral by its type */
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
class LibintIntegralCoreBase : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef typename base_type::scalar_type scalar_type;

  /** Constructor
   *
   * \param system  The molecular system and basis this core deals with
   * \param global  The global libint resources to acquire
   */
  LibintIntegralCoreBase(const LibintSystem& system, const LibintGlobalInit& global)
        : m_system_ptr("LibIntOneElectronIntegralCore", system),
          m_n_bas(static_cast<size_t>(m_system_ptr->basis().nbf())),
          m_global_ptr("LibIntOneElectronIntegralCore", global) {}

  size_t n_rows() const override { return m_n_bas; }
  size_t n_cols() const override { return m_n_bas; }
  size_t n_bas() const { return m_n_bas; }

 protected:
  /** Return the system object */
  const LibintSystem& system() const { return *m_system_ptr; }

 private:
  // Pointer to the molecular system info we use:
  krims::SubscriptionPointer<const LibintSystem> m_system_ptr;

  size_t m_n_bas;  //< Cache for number of basis functions.

  // Pointer to the global libint2 state we need
  krims::SubscriptionPointer<const LibintGlobalInit> m_global_ptr;
};

/** Class for one electron integral cores with libint */
class OneElecIntegralCore final : public LibintIntegralCoreBase {
 public:
  typedef LibintIntegralCoreBase base_type;
  typedef IntegralCoreBase<real_stored_mtx_type> core_base_type;

  /** Constructor
   *
   * \param system  The molecular system and basis this core deals with
   * \param op      The operator libint should compute in this core
   * \param global  The global libint resources to acquire
   */
  OneElecIntegralCore(libint2::Operator op, const LibintSystem& system,
                      const LibintGlobalInit& global)
        : base_type(system, global), m_operator(op) {}

  scalar_type operator()(size_t row, size_t col) const override;

  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A = 1, const scalar_type c_M = 0) const override;

  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

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
  //@{
  /** Call a kernel function for each shell set.
   *
   * Optionally one can supply ranges of rows and columns to work on.
   * Only the shell indices which are within the given ranges are then
   * computed.
   *
   * \param kernel  The kernel function to excecute for each pair of shell indices
   *                This function should take 3 arguments:
   *                        - size_t            Index of the row shell we work on
   *                        - size_t            Index of the col shell we work on
   *                        - const real_type*  Buffer of the integral values.
   * \param row_range   The range of row shell indices to work on (optional)
   * \param col_range   The range of column shell indices to work on (optional)
   **/
  template <typename Kernel>
  void compute_kernel(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
                      Kernel&& kernel) const;

  template <typename Kernel>
  void compute_kernel(Kernel&& kernel) const {
    const auto full_range = krims::range(base_type::system().basis().size());
    compute_kernel(full_range, full_range, std::forward<Kernel>(kernel));
  }
  //@}

 private:
  // Integral operator this core computes.
  libint2::Operator m_operator;
};

class ERICore final : public LibintIntegralCoreBase {
 public:
  typedef LibintIntegralCoreBase base_type;
  typedef IntegralCoreBase<real_stored_mtx_type> core_base_type;
  typedef real_stored_mtx_type stored_matrix_type;
  typedef typename stored_mtx_type::vector_type vector_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  //! Is this exchange or Coulomb operator?
  bool exchange;

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  ERICore(bool exchange_, const LibintSystem& system, const LibintGlobalInit& global)
        : base_type(system, global), exchange(exchange_) {}

  scalar_type operator()(size_t a, size_t b) const override;

  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A = 1, const scalar_type c_M = 0) const override;

  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  std::unique_ptr<core_base_type> clone() const override {
    return std::unique_ptr<core_base_type>(new ERICore(*this));
  }

  IntegralIdentifier id() const override {
    return {IntegralCollection::id,
            (exchange ? IntegralType::exchange : IntegralType::coulomb)};
  }

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap
   **/
  void update(const krims::GenMap& map) override;
};

}  // namespace libint
}  // namespace gaussian
}  // namespace gint
#endif  // GINT_HAVE_LIBINT
