#pragma once
#include "Integral.hh"
#include <complex>
#include <linalgwrap/LazyMatrixExpression.hh>
#include <linalgwrap/type_utils.hh>
#include <string>
#include <type_traits>
#include <vector>

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
  typedef typename stored_matrix_type::scalar_type scalar_type;
  typedef typename linalgwrap::RealTypeOf<scalar_type>::type real_type;
  // end todo

  typedef Integral<stored_matrix_type> integral_matrix_type;

  static_assert(linalgwrap::IsComplexScalar<
                      typename stored_matrix_type::scalar_type>::value &&
                      (otype == COMPLEX_MOLECULAR),
                "The StoredMatrix should only have a complex scalar type if "
                "the OrbitalType is COMPLEX_MOLECULAR");

  /** \name Construct an IntegralLookup object for a particular basis, which
   *  is referred to by the \t basistype_name string */
  IntegralLookup(const std::string& basistype_name);

  /** Return a particular integral, given a name */
  integral_matrix_type operator()(const std::string& integral_name);

  /** Return an integral matrix representing a linear combination of
   *  a list of integrals. */
  auto operator()(const std::vector<scalar_type>& coefficients,
                  const std::vector<std::string>& integral_names)
        // TODO be more clever to figure out the type using decltype syntax here
        -> linalgwrap::LazyMatrixSum<stored_matrix_type>;

private:
  // TODO this is dummy.
  IntegralCoreBase<stored_matrix_type> dummy_integral_core;
};

//
// ---------------------------------------------------
//

template <typename StoredMatrix, OrbitalType otype>
IntegralLookup<StoredMatrix, otype>::IntegralLookup(
      const std::string& basistype_name)
      : dummy_integral_core{} {
  // TODO this dummy stuff;
  std::string local(basistype_name);
}

template <typename StoredMatrix, OrbitalType otype>
typename IntegralLookup<StoredMatrix, otype>::integral_matrix_type
IntegralLookup<StoredMatrix, otype>::operator()(
      const std::string& integral_name) {
  // TODO this dummy stuff;
  std::string local(integral_name);
  return Integral<stored_matrix_type>{dummy_integral_core};
}

template <typename StoredMatrix, OrbitalType otype>
auto IntegralLookup<StoredMatrix, otype>::operator()(
      const std::vector<scalar_type>& coefficients,
      const std::vector<std::string>& integral_names)
      -> linalgwrap::LazyMatrixSum<stored_matrix_type> {
  assert_greater_than(0, coefficients.size());
  assert_size(coefficients.size(), integral_names.size());

  // Iterators to the beginning of vectors:
  auto itcoeff = std::begin(coefficients);
  auto itname = std::begin(integral_names);

  // Note that we have at least one element due to assert_greater_then
  // hence these iterators will always point to an actual element.

  // Lookup the first integral name:
  auto integral0 = (*this)(*itname);

  // TODO get rid of the assumption that the result of the linear combination
  // is a lazy matrix sum?

  // Construct the sum:
  linalgwrap::LazyMatrixSum<stored_matrix_type> sum{std::move(integral0),
                                                    *itcoeff};

  for (; itcoeff != std::end(coefficients); ++itcoeff, ++itname) {
    sum += (*itcoeff) * (*itname);
  }

  return sum;
}
}
