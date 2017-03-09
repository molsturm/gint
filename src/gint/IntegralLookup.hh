#pragma once
#include "Integral.hh"
// TODO It would be nice to get rid of this header since this promotes it to the public
// interface.
#include "IntegralCollectionBase.hh"
#include <linalgwrap/LazyMatrixExpression.hh>
#include <mutex>

namespace gint {

// TODO This orbitaltype interface is awkward and makes coding against this in a
// basis-independent fashion rather annoying.

/** Get Integral objects for a type of basis function */
template <OrbitalType otype>
class IntegralLookup {
 public:
  typedef IntegralCollectionBase<otype> integral_collection_type;
  typedef typename integral_collection_type::stored_mtx_type stored_mtx_type;
  typedef typename stored_mtx_type::scalar_type scalar_type;
  typedef Integral<stored_mtx_type> integral_type;

  static_assert(!krims::IsComplexNumber<scalar_type>::value ||
                      (otype == OrbitalType::COMPLEX_MOLECULAR),
                "The StoredMatrix should only have a complex scalar type if "
                "the OrbitalType is COMPLEX_MOLECULAR");

  /** \name Construct an IntegralLookup object for a particular basis, which
   *  is referred to by the \t basistype_name string
   *
   * ## Some useful control parameters
   *   - "basis_type":    Name of the basis functions to be used.
   *                      (Default: "", which will yield an error)
   *
   *  */
  IntegralLookup(const krims::GenMap& parameters);

  /** Return a particular integral, given an integral type key
   *  (see IntegralType.hh for valid key strings) */
  integral_type lookup_integral(const std::string& integral_type_key) const {
    assert_dbg(m_integral_collection_ptr != nullptr, krims::ExcInternalError());
    return m_integral_collection_ptr->lookup_integral(integral_type_key);
  }

  // TODO later
  // /** Return an integral matrix representing a linear combination of
  //  *  a list of integrals. */
  // integral_type lookup_integral(const std::vector<scalar_type>& coefficients,
  //                               const std::vector<std::string>& integral_type_keys)
  //                               const;

  /** Return the basis type id of the basis type to which this object is i initialised
   *
   * \note Usually this is equivalent to the value supplied as the ``basis_type``
   *       parameter upon construction of this object.
   * */
  const std::string& basis_id() const {
    assert_dbg(m_integral_collection_ptr != nullptr, krims::ExcInternalError());
    return m_integral_collection_ptr->basis_id();
  }

  /** Return the friendly name of the basis type to which this object is initialised */
  const std::string basis_name() const {
    assert_dbg(m_integral_collection_ptr != nullptr, krims::ExcInternalError());
    return m_integral_collection_ptr->basis_name();
  }

  /** Register a basis type with this integral lookup object
   *
   * \param basis_type    The string needed to lookup integrals of this basis later.
   * \param collection_generator   The function which generates an integral collection
   *                               of this type.
   *
   * \note  This function is *not* thread safe
   * */
  static void register_basis_type(std::string basis_type,
                                  collection_generator_type<otype> collection_generator) {
    map_basis_collection_generator[basis_type] = collection_generator;
  }

 private:
  /** The map from the basis id to the generator function for the collection of integrals
   */
  static std::map<std::string, collection_generator_type<otype>>
        map_basis_collection_generator;

  //! Integral collection to use for this lookup class.
  std::unique_ptr<IntegralCollectionBase<otype>> m_integral_collection_ptr;
};

}  // namespace gint
