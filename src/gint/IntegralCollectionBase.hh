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
#include "ERITensor_i.hh"
#include "Integral.hh"
#include "IntegralType.hh"
#include "gint/config.hh"
#include <krims/GenMap.hh>

namespace gint {

DefException1(ExcInvalidIntegralParameters, std::string,
              << "An invalid set of parameters for the integral evaluation was "
                 "encountered:"
              << arg1);

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
 **/
template <typename StoredMatrix>
class IntegralCollectionBase {
 public:
  typedef Integral<StoredMatrix> integral_matrix_type;
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

  /** Get the object representing the repulsion tensor */
  virtual const ERITensor_i<typename StoredMatrix::scalar_type>& eri_tensor() const = 0;

  /** Get the number of basis functions of all integrals returned by this collection. */
  virtual size_t n_bas() const { return lookup_integral(IntegralType::overlap).n_rows(); }

  /** Obtain the id string of the collection / basis type */
  virtual const std::string& basis_id() const = 0;

  /** Obtain the friendly name of the collection / basis type */
  virtual std::string basis_name() const = 0;

  IntegralCollectionBase()                              = default;
  virtual ~IntegralCollectionBase()                     = default;
  IntegralCollectionBase(IntegralCollectionBase&&)      = default;
  IntegralCollectionBase(const IntegralCollectionBase&) = default;
  IntegralCollectionBase& operator=(IntegralCollectionBase&&) = default;
  IntegralCollectionBase& operator=(const IntegralCollectionBase&) = default;
};

//! Type of the generator function, which generates a particular IntegralCollection */
template <typename StoredMatrix>
using collection_generator_type =
      std::function<std::unique_ptr<IntegralCollectionBase<StoredMatrix>>(
            const krims::GenMap& parameters)>;

}  // namespace gint
