#pragma once
#include "Integral.hh"
#include "IntegralCollectionBase.hh"
#include <krims/GenMap.hh>
#include <linalgwrap/LazyMatrixExpression.hh>
#include <type_traits>
#include <vector>

#include "atomic/cs_dummy.hh"
#include "atomic/cs_naive.hh"
#include "atomic/static14.hh"

namespace gint {

static std::map<std::string, create_collection_t<COMPLEX_ATOMIC>*> basis_type_map_ca = {
      {"cs_static14", atomic::static14::IntegralCollection::create},
      {"cs_dummy", atomic::cs_dummy::IntegralCollection::create},
      {"cs_naive", atomic::cs_naive::IntegralCollection::create}};

/** Get Integral objects for a type of basis function */
template <OrbitalType otype>
class IntegralLookup {

 public:
  typedef IntegralCollectionBase<otype> integral_collection_type;
  typedef typename integral_collection_type::stored_mtx_type stored_mtx_type;

  // TODO I feel that these typedefs do not quite belong here ...
  //! The real type underlying the scalar type:
  // TODO make each library know about real type?
  typedef typename stored_mtx_type::scalar_type scalar_type;
  typedef typename stored_mtx_type::real_type real_type;
  // end todo

  static_assert(otype == COMPLEX_ATOMIC, "Currently COMPLEX_ATOMIC is hard coded.");

  // TODO generalise and enable real sturmians
  // typedef typename atomic::static14:IntegralCollection<stored_mtx_type>
  //       integral_collection_type;
  // typedef atomic::cs_dummy::IntegralCollection<stored_mtx_type>
  //      integral_collection_type;
  typedef Integral<stored_mtx_type> integral_type;

  static_assert(!krims::IsComplexNumber<typename stored_mtx_type::scalar_type>::value ||
                      (otype == COMPLEX_MOLECULAR),
                "The StoredMatrix should only have a complex scalar type if "
                "the OrbitalType is COMPLEX_MOLECULAR");

  /** \name Construct an IntegralLookup object for a particular basis, which
   *  is referred to by the \t basistype_name string
   *
   * ## Some useful control parameters
   *   - "basistype_name":    Name of the basis functions to be used.
   *                          (Default: "", which will yield an error)
   *
   *  */
  IntegralLookup(const krims::GenMap& parameters);

  /** Return a particular integral, given an id */
  integral_type operator()(const std::string& integral_id) const;

  /** Return an integral matrix representing a linear combination of
   *  a list of integrals. */
  integral_type operator()(const std::vector<scalar_type>& coefficients,
                           const std::vector<std::string>& integral_names) const;

 private:
  //! Integral collection to use. Currently we use the precomputed static14.
  std::shared_ptr<IntegralCollectionBase<otype>> m_integral_collection;
};

template <OrbitalType otype>
typename IntegralLookup<otype>::integral_type IntegralLookup<otype>::operator()(
      const std::string& integral_id) const {
  return m_integral_collection->lookup_integral(integral_id);
}

template <OrbitalType otype>
typename IntegralLookup<otype>::integral_type IntegralLookup<otype>::operator()(
      const std::vector<scalar_type>& /* coefficients */,
      const std::vector<std::string>& /* integral_names */) const {
  assert_dbg(false, krims::ExcNotImplemented());
  return m_integral_collection("Not implemented");
}

}  // namespace gint
