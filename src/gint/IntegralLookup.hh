#pragma once
#include "Integral.hh"
#include <complex>
#include <linalgwrap/LazyMatrixExpression.hh>
#include <linalgwrap/ParameterMap.hh>
#include <linalgwrap/type_utils.hh>
#include <string>
#include <type_traits>
#include <vector>

#include "atomic/cs_dummy.hh"

namespace gint {

typedef enum {
  REAL_ATOMIC,
  COMPLEX_ATOMIC,
  REAL_MOLECULAR,
  COMPLEX_MOLECULAR
} OrbitalType;

/** Get Integral objects for a type of basis function */
template <typename StoredMatrix, OrbitalType otype>
class IntegralLookup {
public:
  typedef StoredMatrix stored_matrix_type;

  // TODO I feel that these typedefs do not quite belong here ...
  //! The real type underlying the scalar type:
  // TODO make each library know about real type?
  typedef typename stored_matrix_type::scalar_type scalar_type;
  typedef typename linalgwrap::RealTypeOf<scalar_type>::type real_type;
  // end todo

  static_assert(otype == COMPLEX_ATOMIC,
                "Currently only COMPLEX_ATOMIC is hard coded.");
  typedef atomic::cs_dummy::IntegralCollection<stored_matrix_type>
        integral_collection_type;
  typedef Integral<stored_matrix_type> integral_matrix_type;

  static_assert(!linalgwrap::IsComplexScalar<
                      typename stored_matrix_type::scalar_type>::value ||
                      (otype == COMPLEX_MOLECULAR),
                "The StoredMatrix should only have a complex scalar type if "
                "the OrbitalType is COMPLEX_MOLECULAR");

  /** \name Construct an IntegralLookup object for a particular basis, which
   *  is referred to by the \t basistype_name string */
  IntegralLookup(const std::string& basistype_name,
                 const linalgwrap::ParameterMap& parameters =
                       linalgwrap::ParameterMap());

  /** Return a particular integral, given an id */
  integral_matrix_type operator()(const std::string& integral_id);

  /** Return an integral matrix representing a linear combination of
   *  a list of integrals. */
  integral_matrix_type operator()(
        const std::vector<scalar_type>& coefficients,
        const std::vector<std::string>& integral_names);

private:
  //! Integral collection to use. Currently we use the dummy cs_atomic one
  integral_collection_type m_integral_collection;
};

//
// ---------------------------------------------------
//

template <typename StoredMatrix, OrbitalType otype>
IntegralLookup<StoredMatrix, otype>::IntegralLookup(
      const std::string& basistype_name,
      const linalgwrap::ParameterMap& parameters)
      : m_integral_collection{parameters} {
  assert_dbg(basistype_name == "atomic/cs_dummy",
             linalgwrap::ExcNotImplemented());
}

template <typename StoredMatrix, OrbitalType otype>
typename IntegralLookup<StoredMatrix, otype>::integral_matrix_type
IntegralLookup<StoredMatrix, otype>::operator()(
      const std::string& integral_id) {
  return m_integral_collection(integral_id);
}

template <typename StoredMatrix, OrbitalType otype>
typename IntegralLookup<StoredMatrix, otype>::integral_matrix_type
IntegralLookup<StoredMatrix, otype>::operator()(
      const std::vector<scalar_type>& /* coefficients */,
      const std::vector<std::string>& /* integral_names */) {
  assert_dbg(false, linalgwrap::ExcNotImplemented());
  return m_integral_collection("Not implemented");
}

}  // namespace gint
