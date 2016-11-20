#pragma once
#include "Integral.hh"
#include <complex>
#include <krims/ParameterMap.hh>
#include <linalgwrap/LazyMatrixExpression.hh>
#include <string>
#include <type_traits>
#include <vector>

// TODO taken out since it does not work atm
//#include "atomic/cs_dummy.hh"
#include "atomic/static14.hh"

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
  typedef typename stored_matrix_type::real_type real_type;
  // end todo

  static_assert(otype == COMPLEX_ATOMIC,
                "Currently COMPLEX_ATOMIC is hard coded.");

  // TODO generalise and enable real sturmians
  typedef typename atomic::static14::IntegralCollection<stored_matrix_type>
        integral_collection_type;
  // typedef atomic::cs_dummy::IntegralCollection<stored_matrix_type>
  //      integral_collection_type;
  typedef Integral<stored_matrix_type> integral_type;

  static_assert(!krims::IsComplexNumber<
                      typename stored_matrix_type::scalar_type>::value ||
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
  IntegralLookup(const krims::ParameterMap& parameters);

  /** Return a particular integral, given an id */
  integral_type operator()(const std::string& integral_id) const;

  /** Return an integral matrix representing a linear combination of
   *  a list of integrals. */
  integral_type operator()(
        const std::vector<scalar_type>& coefficients,
        const std::vector<std::string>& integral_names) const;

private:
  //! Integral collection to use. Currently we use the precomputed static14.
  integral_collection_type m_integral_collection;
};

//
// ---------------------------------------------------
//

template <typename StoredMatrix, OrbitalType otype>
IntegralLookup<StoredMatrix, otype>::IntegralLookup(
      const krims::ParameterMap& parameters)
      : m_integral_collection{parameters} {
  std::string basis_type =
        parameters.at("basis_type", std::string("No basis_type supplied"));

  assert_dbg(basis_type == integral_collection_type::id,
             krims::ExcNotImplemented());
}

template <typename StoredMatrix, OrbitalType otype>
typename IntegralLookup<StoredMatrix, otype>::integral_type IntegralLookup<
      StoredMatrix, otype>::operator()(const std::string& integral_id) const {
  return m_integral_collection(integral_id);
}

template <typename StoredMatrix, OrbitalType otype>
typename IntegralLookup<StoredMatrix, otype>::integral_type
IntegralLookup<StoredMatrix, otype>::operator()(
      const std::vector<scalar_type>& /* coefficients */,
      const std::vector<std::string>& /* integral_names */) const {
  assert_dbg(false, krims::ExcNotImplemented());
  return m_integral_collection("Not implemented");
}

}  // namespace gint
