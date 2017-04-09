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
#include "IntegralCoreBase.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {
// TODO: Change basis order from n,l,m to m,l,n to make multiplication
// contiguous.

/** Integral core representing the nuclear attraction integral */
class NuclearAttractionIntegralCore final : public IntegralCoreBase {
 public:
  // Compute alpha*A*x + beta*y into y
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** \brief return an element of the matrix \f$ {V_0}_{\mu',\mu} = -Zk/n
   * \delta_{\mu',\mu} \f$ */
  scalar_type operator()(size_t row, size_t col) const override;

  NuclearAttractionIntegralCore(const SturmintSystem& system,
                                const std::string& collection_id)
        : IntegralCoreBase(system, {collection_id, IntegralType::nuclear_attraction}) {}

  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new NuclearAttractionIntegralCore(*this));
  }
};

class KineticIntegralCore final : public IntegralCoreBase {
 public:
  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override;

  KineticIntegralCore(const SturmintSystem& system, const std::string& collection_id)
        : IntegralCoreBase(system, {collection_id, IntegralType::kinetic}) {}

  /** \brief Clone the expression */
  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new KineticIntegralCore(*this));
  }
};

class OverlapIntegralCore final : public IntegralCoreBase {
 public:
  bool has_apply_inverse() const final override { return true; }

  void apply(const const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  void apply_inverse(const const_multivector_type& x, multivector_type& y,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_t row, size_t col) const override;

  OverlapIntegralCore(const SturmintSystem& system, const std::string& collection_id)
        : IntegralCoreBase(system, {collection_id, IntegralType::overlap}) {}

  /** \brief Clone the expression */
  std::unique_ptr<base_core_type> clone() const override {
    return std::unique_ptr<base_core_type>(new OverlapIntegralCore(*this));
  }
};

//
// ------------------------------------------------------------------------------
//

// Preamble for apply calls. Not that all modes are the same since we are symmetric
// and real, so we do not switch over mode
#define APPLY_PRECONDITIONS()                                              \
  {                                                                        \
    assert_finite(c_A);                                                    \
    assert_finite(c_y);                                                    \
    assert_size(x.n_cols(), y.n_cols());                                   \
    assert_size(x.n_rows(), n_cols());                                     \
    assert_size(y.n_rows(), n_rows());                                     \
    assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans); \
  }

inline void NuclearAttractionIntegralCore::apply(const const_multivector_type& x,
                                                 multivector_type& y,
                                                 const linalgwrap::Transposed mode,
                                                 const scalar_type c_A,
                                                 const scalar_type c_y) const {
  APPLY_PRECONDITIONS();

  const scalar_type* x_ptr = x.data().memptr();
  scalar_type* y_ptr = const_cast<scalar_type*>(y.data().memptr());
  const int x_cols = static_cast<int>(x.n_cols());

  const auto Z = system().Z;
  const auto k = system().k;

  using sturmint::atomic::cs::apply_to_full_vectors::nuclear_attraction;
  nuclear_attraction<real_type>(system().basis, x_ptr, y_ptr, Z * k * c_A, c_y, x_cols);
}

inline scalar_type NuclearAttractionIntegralCore::operator()(size_t row,
                                                             size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());
  const int n = system().basis[row].n;

  return row == col ? -system().Z * system().k / n : 0;
}

//

inline void KineticIntegralCore::apply(const const_multivector_type& x,
                                       multivector_type& y,
                                       const linalgwrap::Transposed mode,
                                       const scalar_type c_A,
                                       const scalar_type c_y) const {
  APPLY_PRECONDITIONS();

  const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
  real_type* y_ptr = const_cast<real_type*>(y.data().memptr());
  const int x_cols = static_cast<int>(x.n_cols());

  const auto k = system().k;
  const scalar_type c_kkA = static_cast<scalar_type>(-0.5L * k * k * c_A);

  using sturmint::atomic::cs::apply_to_full_vectors::overlap;
  overlap<scalar_type>(system().basis, x_ptr, y_ptr, c_kkA, c_y, x_cols);
  y += c_A * k * k * x;  // kinetic(x) = k^2*x-1/2 overlap(x)
}

inline scalar_type KineticIntegralCore::operator()(size_t row, size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  const nlm_type mui = system().basis[row];
  const nlm_type muj = system().basis[col];
  const auto k = system().k;
  return static_cast<scalar_type>(k * k * sturmint::atomic::cs::kinetic(mui, muj));
}

//

inline void OverlapIntegralCore::apply(const const_multivector_type& x,
                                       multivector_type& y,
                                       const linalgwrap::Transposed mode,
                                       const scalar_type c_A,
                                       const scalar_type c_y) const {
  APPLY_PRECONDITIONS();

  const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
  real_type* y_ptr = const_cast<real_type*>(y.data().memptr());
  const int x_cols = static_cast<int>(x.n_cols());

  using sturmint::atomic::cs::apply_to_full_vectors::overlap;
  overlap(system().basis, x_ptr, y_ptr, c_A, c_y, x_cols);
}

inline void OverlapIntegralCore::apply_inverse(const const_multivector_type& x,
                                               multivector_type& y,
                                               const linalgwrap::Transposed mode,
                                               const scalar_type c_A,
                                               const scalar_type c_y) const {
  APPLY_PRECONDITIONS();

  const real_type* x_ptr = const_cast<const real_type*>(x.data().memptr());
  real_type* y_ptr = const_cast<real_type*>(y.data().memptr());

  using sturmint::atomic::cs::apply_to_full_vectors::overlap_inverse;
  overlap_inverse(system().basis, x_ptr, y_ptr, c_A, c_y, static_cast<int>(x.n_cols()));
}

inline scalar_type OverlapIntegralCore::operator()(size_t row, size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  const nlm_type mui = system().basis[row];
  const nlm_type muj = system().basis[col];
  return static_cast<scalar_type>(sturmint::atomic::cs::overlap(mui, muj));
}

#undef APPLY_PRECONDITIONS

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
