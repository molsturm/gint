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
#ifdef GINT_HAVE_STATIC_INTEGRALS

#include "Static14Data.hh"
#include "gint/CoefficientContainer.hh"
#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/config.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_static14 {

// In this namespace all things are real:
using namespace gint::real_valued;

void apply_stored_matrix(const stored_matrix_type& A, const const_multivector_type& x,
                         multivector_type& y,
                         const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                         const scalar_type c_A = 1, const scalar_type c_y = 0);

class ERITensor final : public ERITensor_i<scalar_type> {
 public:
  const real_type k;

  ERITensor(real_type k) : k(k) {}

  size_t n_bas() const override { return Static14Data::nbas; }

 protected:
  using ERITensor_i<scalar_type>::compute_kernel;
  void compute_kernel(const std::array<krims::Range<size_t>, 4>& block,
                      kernel_type kernel) const override;
};

/** This integral collection contains precomputed data for an atomic sturmian
 * basis with n_max = 3 and l_max=2, which yields 14 basis functions.
 *
 * The exponent k and the nuclear charge Z can still be chosen freely.
 */
class IntegralCollection final : public IntegralCollectionBase<stored_matrix_type> {
 public:
  typedef IntegralCollectionBase<stored_matrix_type> base_type;

  const static std::string id;
  const real_type k_exponent, Z_charge;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - structure (gint::Structure): The structure of the system (i.e. an atom)
   */
  IntegralCollection(const krims::GenMap& parameters);

  using base_type::lookup_integral;
  integral_matrix_type lookup_integral(IntegralType type) const override;

  const ERITensor_i<scalar_type>& eri_tensor() const override { return m_eri_tensor; }

  /** Obtain the id string of the collection / basis type */
  const std::string& basis_id() const override { return id; }

  /** Obtain the friendly name of the collection / basis type */
  std::string basis_name() const override {
    return "Fully precomputed integral data for 14 atomic coulomb sturmians";
  }

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::unique_ptr<base_type> create(const krims::GenMap& parameters) {
    return krims::make_unique<IntegralCollection>(parameters);
  }

 private:
  ERITensor m_eri_tensor;
};

//
// Integral cores
//
class OneElecIntegralCore final : public gint::IntegralCoreBase<stored_matrix_type> {
 public:
  typedef gint::IntegralCoreBase<stored_matrix_type> base_type;
  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return Static14Data::nbas; }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return Static14Data::nbas; }

  /** \brief return an element of the operator matrix */
  scalar_type operator()(size_t row, size_t col) const override {
    return m_fac * (*m_mat_ptr)(row, col);
  }

  bool has_transpose_operation_mode() const override {
    return m_mat_ptr->has_transpose_operation_mode();
  }

  bool has_apply_inverse() const override { return m_inv_mat_ptr != nullptr; }

  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override {
    apply_stored_matrix(*m_mat_ptr, x, y, mode, m_fac * c_A, c_y);
  }

  void apply_inverse(const const_multivector_type& x, multivector_type& y,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A             = 1,
                     const scalar_type c_y             = 0) const override {
    assert_throw(m_inv_mat_ptr,
                 krims::ExcDisabled("The apply_inverse function is in general "
                                    "very expensive and is only implemented in "
                                    "some cases. Use the function "
                                    "has_apply_inverse() to check when."));
    apply_stored_matrix(*m_inv_mat_ptr, x, y, mode, c_A / m_fac, c_y);
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   */
  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_this          = 1,
                     const scalar_type c_M             = 0) const override {
    m_mat_ptr->extract_block(M, start_row, start_col, mode, m_fac * c_this, c_M);
  }

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new OneElecIntegralCore(*this));
  }

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override { return {IntegralCollection::id, m_type}; }

  OneElecIntegralCore(IntegralType type, const scalar_type Z, const scalar_type k);

 private:
  //! Factor to implicitly multiply the static data with
  scalar_type m_fac;

  //! Inner matrix
  krims::SubscriptionPointer<const stored_matrix_type> m_mat_ptr;

  //! Inverse matrix
  krims::SubscriptionPointer<const stored_matrix_type> m_inv_mat_ptr;

  //! Integral type
  IntegralType m_type;
};

class ERICore final : public gint::IntegralCoreBase<stored_matrix_type>,
                      public CoefficientContainer<stored_matrix_type> {
 public:
  typedef gint::IntegralCoreBase<stored_matrix_type> base_type;

  //! The exponent
  const real_type k;

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return Static14Data::nbas; }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return Static14Data::nbas; }

  /** \brief return an element of the matrix \f$ {J}_{\mu',\mu} \f$ or \f$
   * {k}_{\mu',\mu} \f$. */
  scalar_type operator()(size_t a, size_t b) const override;

  bool has_transpose_operation_mode() const override { return true; }

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new ERICore(*this));
  }

  void update(const krims::GenMap& map) override {
    CoefficientContainer<stored_matrix_type>::update(map);
    assert_dbg(coeff_bo_ptr == nullptr || coeff_bo().n_vectors() == 0 ||
                     coeff_bo_ptr->n_elem() == Static14Data::nbas,
               krims::ExcSizeMismatch(coeff_bo_ptr->n_elem(), Static14Data::nbas));
  }

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override { return {IntegralCollection::id, type}; }

  ERICore(IntegralType type, real_type k) : k(k), type(type) {}

 private:
  IntegralType type;
};

}  // namespace cs_static14
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint

#endif  // GINT_HAVE_STATIC_INTEGRALS
