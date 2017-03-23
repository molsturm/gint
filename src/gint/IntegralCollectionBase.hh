#pragma once
#include "Integral.hh"
#include "IntegralType.hh"
#include "gint/config.hh"
#include <krims/GenMap.hh>

namespace gint {

enum class OrbitalType {
  /** Real atomic orbitals */
  REAL_ATOMIC,

  /** Complex atomic orbitals */
  COMPLEX_ATOMIC,

  /** Real molecular orbitals */
  REAL_MOLECULAR,

  /** Complex molecular orbitals */
  COMPLEX_MOLECULAR
};

/** Base class for managing collections of integrals
 *
 * Next to the interface enforced by the pure virtual methods in this class
 * each derived class should contain the following static members:
 *    - A function ``std::unique_ptr<IntegralCollectionBase<otype>> create(
 *      const krims::GenMap&)``
 *      which creates an object of said collection using the provided parameters.
 *    - A ``const std::string id`` which is used as the id of the basis /
 *      integral collection. A user of the IntegralLookup library will need these
 *      ids as keys to request the integrals in the said collection
 *      (This should be returned by basis_id)
 * */
template <OrbitalType otype>
class IntegralCollectionBase {
 public:
  constexpr static OrbitalType orbital_type = otype;
  typedef typename std::conditional<otype == OrbitalType::COMPLEX_MOLECULAR,
                                    detail::complex_stored_mtx_type,
                                    detail::real_stored_mtx_type>::type stored_mtx_type;
  typedef Integral<stored_mtx_type> integral_matrix_type;

  /** Lookup an integral in this collection by its integral type key
   *
   * \throws ExcInvalidIntegralTypeKey if the integral type is not valid.
   * */
  virtual integral_matrix_type lookup_integral(
        const std::string& integral_type_key) const {
    return lookup_integral(to_integral_type(integral_type_key));
  }

  /** Lookup an integral in this collection by its integral type */
  virtual integral_matrix_type lookup_integral(IntegralType type) const = 0;

  /** Obtain the id string of the collection / basis type */
  virtual const std::string& basis_id() const = 0;

  /** Obtain the friendly name of the collection / basis type */
  virtual std::string basis_name() const = 0;

  virtual ~IntegralCollectionBase() = default;
  IntegralCollectionBase& operator=(const IntegralCollectionBase&) = default;
  IntegralCollectionBase& operator=(IntegralCollectionBase&&) = default;
  IntegralCollectionBase(const IntegralCollectionBase&) = default;
  IntegralCollectionBase(IntegralCollectionBase&&) = default;
  IntegralCollectionBase() = default;
};

//! Type of the generator function, which generates a particular IntegralCollection */
template <OrbitalType otype>
using collection_generator_type =
      std::function<std::unique_ptr<IntegralCollectionBase<otype>>(
            const krims::GenMap& parameters)>;

}  // namespace gint
