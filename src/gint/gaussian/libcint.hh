//
// Copyright (C) 2017 by the gint authors
//
// This file is part of gint.
//
// gint is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gint is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gint. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "gint/config.hh"
#ifdef GINT_HAVE_LIBCINT

#include <krims/SubscriptionPointer.hh>

#include "Basis.hh"
#include "gint/CoefficientContainer.hh"
#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/Structure.hh"

namespace gint {
namespace gaussian {
namespace libcint {

// In this namespace all things are real:
using namespace real_valued;

// The integer type used by libcint
typedef int int_type;

/** Object representing a molecular structure and basis in a way usable to libcint.
 *
 * In a way this is the system (physical and model) we want to compute integrals for.
 *
 * \note  Modifying this object after construction is not supported since many of
 *        the returned properties are cached values which are never invalidated.
 **/
class System : public krims::Subscribable {
 public:
  System() : shell_data{}, atom_data{}, m_n_bas(0), m_n_shells(0), m_n_atoms(0) {}
  System(krims::SubscriptionPointer<const Structure> structure_ptr, Basis basis);

  static constexpr size_t env_size = 10000;

  /** Environment which stores the coordinates, atom positions and contraction
   *  coefficients and exponents **/
  std::array<double, env_size> env{};

  /** Container holding the basis information more precisely the
   *  information about the shells of the basis.
   *
   *  \note This thing is called "bas" inside libint.
   */
  std::vector<int_type> shell_data;

  //! Container holding the atom information
  std::vector<int_type> atom_data;

  /** Access to the total number of basis functions */
  size_t n_bas() const { return m_n_bas; }

  /** Return the number of shells
   *
   * \note This is *not* the size of the shell_data array!
   * */
  size_t n_shells() const { return m_n_shells; }

  /** Return the number of atoms
   *
   * \note This is *not* the size of the atom_data array!
   */
  size_t n_atoms() const { return m_n_atoms; }

 private:
  static_assert(env_size < std::numeric_limits<int>::max() - 5,
                "The number of elements in env needs to be representable by an int.");

  //! The index offset to the next free slot in the env array
  size_t env_back_off() const { return static_cast<size_t>(env_back - env.begin()); }

  /** The number of basis functions */
  size_t m_n_bas;

  /** The number of shells for which information is stored in the shells array */
  size_t m_n_shells;

  /** The number of atoms for which information is stored in the atoms array */
  size_t m_n_atoms;

  //! Pointer to the next free slot in the env array.
  typename std::array<double, env_size>::iterator env_back = env.begin();
};

/** Repulsion integral tensor object for libcint. */
class ERITensor final : public ERITensor_i<scalar_type> {
 public:
  /** Constructor
   * \param system  The molecular system and basis this core deals with
   */
  ERITensor(const System& system) : m_system_ptr("Libcint::ERITensor", system) {}

  /** Number of basis functions in this bases */
  size_t n_bas() const override { return m_system_ptr->n_bas(); }

 protected:
  using ERITensor_i::compute_kernel;

  void compute_kernel(const std::array<krims::Range<size_t>, 4>& block,
                      kernel_type kernel) const override;

 private:
  // Pointer to the molecular system info we use:
  krims::SubscriptionPointer<const System> m_system_ptr;
};

/** IntegralCollection for the libcint integral objects */
class IntegralCollection final : public IntegralCollectionBase<stored_matrix_type> {
 public:
  typedef IntegralCollectionBase<stored_matrix_type> base_type;

  const static std::string id;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - basis_set  (std::string)
   *     Either the name of a Gaussian basis set to use or the
   *     path to a basis set file to read the basis set from.
   *   - basis     (gint::gaussian::Basis)  The basis to use to model the
   *                                        molecule.
   *   - structure (gint::Structure)  The structure of the molecule.
   */
  explicit IntegralCollection(const krims::GenMap& parameters);

  using base_type::lookup_integral;
  integral_matrix_type lookup_integral(IntegralType type) const override;

  const ERITensor_i<scalar_type>& eri_tensor() const override { return m_eri_tensor; }

  /** Obtain the id string of the collection / basis type */
  const std::string& basis_id() const override { return id; }

  /** Return the number of basis functions */
  size_t n_bas() const override { return m_system.n_bas(); }

  /** Obtain the friendly name of the collection / basis type */
  std::string basis_name() const override { return "Gaussian integrals from libcint"; }

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::unique_ptr<base_type> create(const krims::GenMap& parameters) {
    return krims::make_unique<IntegralCollection>(parameters);
  }

 private:
  const System m_system;
  ERITensor m_eri_tensor;
};

//
// Integral Cores
//

/** Base class for libint integral cores */
class LibCintIntegralCoreBase : public IntegralCoreBase<stored_matrix_type> {
 public:
  typedef IntegralCoreBase<stored_matrix_type> base_type;
  typedef typename base_type::scalar_type scalar_type;

  /** Constructor
   *
   * \param system  The molecular system and basis this core deals with
   * \param global  The global libint resources to acquire
   */
  LibCintIntegralCoreBase(const System& system)
        : m_system_ptr("LibIntOneElectronIntegralCore", system) {}

  size_t n_rows() const override { return m_system_ptr->n_bas(); }
  size_t n_cols() const override { return m_system_ptr->n_bas(); }
  size_t n_bas() const { return m_system_ptr->n_bas(); }

  scalar_type operator()(size_t row, size_t col) const override;

  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const lazyten::Transposed mode = lazyten::Transposed::None,
                     const scalar_type c_A = 1, const scalar_type c_M = 0) const override;

  void apply(const const_multivector_type& x, multivector_type& y,
             const lazyten::Transposed mode = lazyten::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

 protected:
  /** Return the system object */
  const System& system() const { return *m_system_ptr; }

  /** Type of a computational kernel for the compute function below.
   *
   * For each pair of shell indices (e.g. 1s, 2s, 2p, ...) the kernel function
   * is called once.
   */
  typedef std::function<void(size_t, size_t, const scalar_type*)> kernel_type;

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
    assert_internal(n_rows() == n_cols());
    const auto full_range = krims::range(n_rows());
    compute(full_range, full_range, std::forward<kernel_type>(kernel));
  }
  //@}

 private:
  // Pointer to the molecular system info we use:
  krims::SubscriptionPointer<const System> m_system_ptr;
};

/** Class for one electron integral cores with libint */
class OneElecIntegralCore final : public LibCintIntegralCoreBase {
 public:
  typedef LibCintIntegralCoreBase base_type;
  typedef IntegralCoreBase<stored_matrix_type> core_base_type;

  /** Constructor
   *
   * \param system  The molecular system and basis this core deals with
   * \param op      The operator libint should compute in this core
   * \param global  The global libint resources to acquire
   */
  OneElecIntegralCore(IntegralType type, const System& system);

  std::unique_ptr<core_base_type> clone() const override {
    return std::unique_ptr<core_base_type>(new OneElecIntegralCore(*this));
  }

  IntegralIdentifier id() const override { return {IntegralCollection::id, m_type}; }

 protected:
  virtual void compute(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
                       typename base_type::kernel_type&& kernel) const override;

 private:
  int_type (*m_cint1e_kernel)(double*, const int_type*, const int_type*, const int_type,
                              const int_type*, const int_type, const double*);

  // Integral operator this core computes.
  const IntegralType m_type;
};

class ERICore final : public LibCintIntegralCoreBase,
                      public CoefficientContainer<stored_matrix_type> {
 public:
  typedef LibCintIntegralCoreBase base_type;
  typedef IntegralCoreBase<stored_matrix_type> core_base_type;

  ERICore(IntegralType type, const System& system) : base_type(system), m_type(type) {
    assert_internal(type == IntegralType::exchange || type == IntegralType::coulomb);
  }

  std::unique_ptr<core_base_type> clone() const override {
    return std::unique_ptr<core_base_type>(new ERICore(*this));
  }

  IntegralIdentifier id() const override { return {IntegralCollection::id, m_type}; }

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap
   **/
  void update(const krims::GenMap& map) override {
    CoefficientContainer<stored_matrix_type>::update(map);
    assert_dbg(coeff_bo_ptr == nullptr || coeff_bo().n_vectors() == 0 ||
                     coeff_bo_ptr->n_elem() == n_bas(),
               krims::ExcSizeMismatch(coeff_bo_ptr->n_elem(), n_bas()));
  }

 protected:
  void compute(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
               kernel_type&& kernel) const override;

  IntegralType m_type;
};

}  // namespace libcint
}  // namespace gaussian
}  // namespace gint

#endif  // GINT_HAVE_LIBCINT
