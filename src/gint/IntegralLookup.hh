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
#include "Integral.hh"
#include "IntegralCollectionBase.hh"

namespace gint {

/** Get Integral objects for a type of basis function */
template <typename StoredMatrix>
class IntegralLookup {
 public:
  typedef StoredMatrix stored_matrix_type;
  typedef Integral<stored_matrix_type> integral_type;

  /** \name Construct an IntegralLookup object for a particular basis, which
   *  is referred to by the \t basistype_name string
   *
   * ## Some useful control parameters
   *   - "basis_type":    Name of the basis functions to be used.
   *                      (Default: "", which will yield an error)
   *   - "orbital_type":  The type of orbitals to be used.
   *                      (Default: REAL_MOLECULAR )
   **/
  IntegralLookup(const krims::GenMap& parameters);

  /** Return a particular integral, given an integral type key
   *  (see IntegralType.hh for valid key strings) */
  integral_type lookup_integral(const std::string& integral_type_key) const {
    assert_internal(m_integral_collection_ptr != nullptr);
    return m_integral_collection_ptr->lookup_integral(integral_type_key);
  }

  // TODO later
  // /** Return an integral matrix representing a linear combination of
  //  *  a list of integrals. */
  // integral_type lookup_integral(const std::vector<scalar_type>& coefficients,
  //                               const std::vector<std::string>& integral_type_keys)
  //                               const;

  const ERITensor_i<typename StoredMatrix::scalar_type>& eri_tensor() const {
    assert_internal(m_integral_collection_ptr != nullptr);
    return m_integral_collection_ptr->eri_tensor();
  }

  /* Return the number of basis functions which are used to discretise the
   * operators in this IntegralLookup. In other words this is the number
   * of rows and columns of (most of) the integrals like the overlap and so on.
   */
  size_t n_bas() const {
    assert_internal(m_integral_collection_ptr != nullptr);
    return m_integral_collection_ptr->n_bas();
  }

  /** Return the basis type id of the basis type to which this object is i initialised
   *
   * \note Usually this is equivalent to the value supplied as the ``basis_type``
   *       parameter upon construction of this object.
   * */
  const std::string& basis_id() const {
    assert_internal(m_integral_collection_ptr != nullptr);
    return m_integral_collection_ptr->basis_id();
  }

  /** Return the friendly name of the basis type to which this object is initialised */
  const std::string basis_name() const {
    assert_internal(m_integral_collection_ptr != nullptr);
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
  static void register_basis_type(
        std::string basis_type,
        collection_generator_type<StoredMatrix> collection_generator) {
    map_basis_collection_generator[basis_type] = collection_generator;
  }

  /** Return the list of basis types which are registered at the very moment
   *
   * \note  This function is not thread safe
   */
  static std::vector<std::string> available_basis_types();
  // TODO The above function should also exist in a way which does not distinguish
  //      between real and complex.

 private:
  /** The map from the basis id to the generator function for the collection
   *  of integrals */
  static std::map<std::string, collection_generator_type<StoredMatrix>>
        map_basis_collection_generator;

  //! Integral collection to use for this lookup class.
  std::unique_ptr<IntegralCollectionBase<StoredMatrix>> m_integral_collection_ptr;
};

}  // namespace gint
