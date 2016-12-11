#pragma once
#ifdef GINT_STATIC_INTEGRALS

#include "Static14Data.hh"
#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/config.hh"
#include <krims/ParameterMap.hh>

namespace gint {
namespace atomic {
namespace static14 {

class OverlapIntegralCore;
class NuclearAttractionIntegralCore;
class KineticIntegralCore;
class ERICore;

#include "gint/real_config.hh"

void apply_stored_matrix(const real_stored_mtx_type& A, const_multivector_type& x,
                         multivector_type& y,
                         const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                         const scalar_type c_A = 1, const scalar_type c_y = 0);

/** This integral collection contains precomputed data for an atomic sturmian
 * basis with n_max = 3 and l_max=2, which yields 14 basis functions.
 *
 * The exponent k and the nuclear charge Z can still be chosen freely.
 */
class IntegralCollection : public IntegralCollectionBase<COMPLEX_ATOMIC> {
 public:
  typedef IntegralCollectionBase<COMPLEX_ATOMIC> base_type;

  const static std::string id, name;
  const real_type k_exponent, Z_charge;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - Z_charge (double): The nuclear change of the system
   */
  IntegralCollection(const krims::ParameterMap& parameters);

  /** Lookup an integral by its identifier string */
  integral_matrix_type lookup_integral(const std::string& integral_id) const override;

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::shared_ptr<base_type> create(const krims::ParameterMap& parameters) {
    return std::make_shared<IntegralCollection>(parameters);
  }
};

class NuclearAttractionIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef real_stored_mtx_type stored_mtx_type;
  typedef IntegralCoreBase<real_stored_mtx_type> base_type;
  typedef real_type scalar_type;

  const real_type k, Z;

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief return an element of the matrix \f$ {V_0}_{\mu',\mu} = -Zk/n
   * \delta_{\mu',\mu} \f$ */
  scalar_type operator()(size_t row, size_t col) const override {
    return -k * Z * detail::Static14Data<stored_mtx_type>::v0_bb_base(row, col);
  }

  bool has_transpose_operation_mode() const override {
    return detail::Static14Data<stored_mtx_type>::v0_bb_base
          .has_transpose_operation_mode();
  }

  // c_A * A * x + c_y * y
  void apply(const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override {
    apply_stored_matrix(detail::Static14Data<stored_mtx_type>::v0_bb_base, x, y, mode,
                        -k * Z * c_A, c_y);
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   */
  virtual void extract_block(
        stored_matrix_type& M, const size_t start_row, const size_t start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M = linalgwrap::Constants<scalar_type>::zero) const override {
    detail::Static14Data<stored_matrix_type>::v0_bb_base.extract_block(
          M, start_row, start_col, mode, -k * Z * c_this, c_M);
  }

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new NuclearAttractionIntegralCore(*this));
  }

  // TODO use static keys to generate values
  /** \brief Get the identifier of the integral */
  std::string id() const override { return "static/14/nuclear_attraction"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Nuclear attraction operator"; }

  NuclearAttractionIntegralCore(real_type k, real_type Z) : k(k), Z(Z) {}
};

//
// Integral cores
//
class OverlapIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef real_stored_mtx_type stored_mtx_type;
  typedef IntegralCoreBase<stored_mtx_type> base_type;
  typedef real_type scalar_type;

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief return an element of the matrix \f$ {S}_{\mu',\mu} */
  scalar_type operator()(size_t row, size_t col) const override {
    return detail::Static14Data<stored_mtx_type>::s_bb(row, col);
  }

  bool has_transpose_operation_mode() const override {
    return detail::Static14Data<stored_mtx_type>::s_bb.has_transpose_operation_mode();
  }

  bool has_apply_inverse() const override { return true; }

  void apply(const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override {
    apply_stored_matrix(detail::Static14Data<stored_mtx_type>::s_bb, x, y, mode, c_A,
                        c_y);
  }

  void apply_inverse(const_multivector_type& x, multivector_type& y,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_A = 1,
                     const scalar_type c_y = 0) const override {
    apply_stored_matrix(detail::Static14Data<stored_mtx_type>::sinv_bb, x, y, mode, c_A,
                        c_y);
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.  */
  virtual void extract_block(
        stored_matrix_type& M, const size_t start_row, const size_t start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M = linalgwrap::Constants<scalar_type>::zero) const override {
    detail::Static14Data<stored_matrix_type>::s_bb.extract_block(M, start_row, start_col,
                                                                 mode, c_this, c_M);
  }

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new OverlapIntegralCore(*this));
  }

  // TODO use static keys to generate values
  /** \brief Get the identifier of the integral */
  std::string id() const override { return "static/14/overlap"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Overlap operator"; }
};

class KineticIntegralCore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef real_stored_mtx_type stored_mtx_type;
  typedef IntegralCoreBase<stored_mtx_type> base_type;
  typedef real_type scalar_type;

  const real_type k;

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief return an element of the matrix \f$ {T}_{\mu',\mu} */
  scalar_type operator()(size_t row, size_t col) const override {
    return k * k * detail::Static14Data<stored_mtx_type>::t_bb_base(row, col);
  }

  bool has_transpose_operation_mode() const override {
    return detail::Static14Data<stored_mtx_type>::t_bb_base
          .has_transpose_operation_mode();
  }

  void apply(const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override {
    apply_stored_matrix(detail::Static14Data<stored_mtx_type>::t_bb_base, x, y, mode,
                        k * k * c_A, c_y);
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.  */
  virtual void extract_block(
        stored_matrix_type& M, const size_t start_row, const size_t start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M = linalgwrap::Constants<scalar_type>::zero) const override {
    detail::Static14Data<stored_matrix_type>::t_bb_base.extract_block(
          M, start_row, start_col, mode, k * k * c_this, c_M);
  }

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new KineticIntegralCore(*this));
  }

  // TODO use static keys to generate values
  /** \brief Get the identifier of the integral */
  std::string id() const override { return "static/14/kinetic"; }

  /** \brief Get the friendly name of the integral */
  std::string name() const override { return "Kinetic energy operator"; }

  KineticIntegralCore(real_type k) : k(k) {}
};

class ERICore : public IntegralCoreBase<real_stored_mtx_type> {
 public:
  typedef real_stored_mtx_type stored_mtx_type;
  typedef IntegralCoreBase<stored_mtx_type> base_type;
  typedef typename stored_mtx_type::vector_type vector_type;
  typedef real_type scalar_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  //! Is this an exchange or Coulomb operator?
  const bool exchange;

  //! The exponent
  const real_type k;

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  /** \brief Number of rows of the matrix */
  size_t n_rows() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief Number of columns of the matrix  */
  size_t n_cols() const override { return detail::Static14Data<stored_mtx_type>::nbas; }

  /** \brief return an element of the matrix \f$ {J}_{\mu',\mu} \f$ or \f$
   * {k}_{\mu',\mu} \f$. */
  scalar_type operator()(size_t a, size_t b) const override {
    assert_greater(a, n_rows());
    assert_greater(b, n_cols());

    // This routine computes
    //
    // J_{ab} = \sum_{cd} P_{cd} < ac | bd >
    //        = \sum_{cd} \sum_{k \in occ} (C^{k})_c (C^{k})_d < ac | bd >
    // where a and b are the same centre, so are c and d
    //
    // or alternatively
    //
    // K_{ab} = \sum_{cd} P_{cd} < ab | cd >
    // where a and c are the same centre, so are b and d

    // Repulsion integrals indexed in shell-pairs
    const auto& i_bbbb = detail::Static14Data<stored_mtx_type>::i_bbbb_base;
    const size_t nbas = detail::Static14Data<stored_mtx_type>::nbas;

    assert_dbg(coefficients_occupied_ptr != nullptr, krims::ExcInvalidPointer());

    // Density matrix expression
    auto density_bb =
          outer_prod_sum(*coefficients_occupied_ptr, *coefficients_occupied_ptr);
    assert_dbg(density_bb.n_rows() == nbas, krims::ExcInternalError());
    assert_dbg(density_bb.n_cols() == nbas, krims::ExcInternalError());

    // Shell pair index for basis functions a and b:
    const size_t ab_pair = a * nbas + b;

    // Sum accumulator variable for this exchange or
    // coulomb  matrix element
    scalar_type mat_ab{0};

    // Double loop over basis functions c and d:
    for (size_t c = 0; c < nbas; ++c) {
      // Shell pair index for basis functions a and c:
      const size_t ac_pair = a * nbas + c;

      for (size_t d = 0; d < nbas; ++d) {
        // Shell pair index for basis functions c and d:
        // or b and d:
        const size_t cd_pair = c * nbas + d;
        const size_t db_pair = d * nbas + b;

        // Perform contraction:
        const scalar_type i_elem =
              exchange ? i_bbbb(ac_pair, db_pair) : i_bbbb(ab_pair, cd_pair);
        mat_ab += k * i_elem * density_bb(c, d);
      }  // d
    }    // c

    return mat_ab;
  }

  bool has_transpose_operation_mode() const override { return true; }

  void apply(
        const_multivector_type& x, multivector_type& y,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_A = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_y = linalgwrap::Constants<scalar_type>::zero) const override {
    using namespace linalgwrap;
    const size_t nbas = detail::Static14Data<stored_mtx_type>::nbas;

    assert_finite(c_A);
    assert_finite(c_y);
    assert_size(x.n_cols(), y.n_cols());
    assert_size(x.n_rows(), nbas);
    assert_size(y.n_rows(), nbas);
    assert_sufficiently_tested(mode != Transposed::ConjTrans);
    // All modes are same case since we are symmetric and real, so no
    // switching over mode.

    // Scale the current values of y or set them to zero
    // (if c_y == 0): We are now done with c_y and do not
    // need to worry about it any more in this function
    for (size_t i = 0; i < y.n_rows(); i++)
      for (size_t j = 0; j < y.n_cols(); j++) y(i, j) = (c_y != 0 ? c_y * y(i, j) : 0);

    // if c_this == 0 we are done
    if (c_A == Constants<scalar_type>::zero) return;

    for (size_t veci = 0; veci < x.n_cols(); ++veci) {
      for (size_t row = 0; row < nbas; ++row) {
        scalar_type sum = 0;
        for (size_t k = 0; k < nbas; ++k) sum += (*this)(row, k) * x(k, veci);

        y(row, veci) += c_A * sum;
      }  // row
    }    // veci
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.  */
  virtual void extract_block(
        stored_matrix_type& M, const size_t start_row, const size_t start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M = linalgwrap::Constants<scalar_type>::zero) const override {
    using namespace linalgwrap;
#ifdef DEBUG
    const size_t nbas = detail::Static14Data<stored_matrix_type>::nbas;
#endif

    assert_finite(c_this);
    assert_finite(c_M);
    assert_greater_equal(start_row + M.n_rows(), nbas);
    assert_greater_equal(start_col + M.n_cols(), nbas);
    assert_sufficiently_tested(mode != Transposed::ConjTrans);

    // For empty matrices there is nothing to do
    if (M.n_rows() == 0 || M.n_cols() == 0) return;

    // Scale the current values of M or set them to zero
    // (if c_M == 0): We are now done with c_M and do not
    // need to worry about it any more in this function
    // TODO we are kind of calling an internal function here
    linalgwrap::detail::scale_or_set(M, c_M);

    if (c_this == Constants<scalar_type>::zero) return;

    for (size_t row = 0; row < M.n_rows(); ++row) {
      for (size_t col = 0; col < M.n_cols(); ++col) {
        switch (mode) {
          case Transposed::None:
          case Transposed::Trans:
            // No difference since symmetric
            M(row, col) += c_this * (*this)(start_row + row, start_col + col);
            break;
          case Transposed::ConjTrans:
            // A variant of std::conj, which does not return a complex
            // data type if scalar is real only.
            // TODO we are kind of calling an internal function here
            linalgwrap::detail::ConjFctr mconj;
            M(row, col) += c_this * mconj((*this)(start_col + col, start_row + row));
            break;
        }  // mode
      }    // col
    }      // row
  }

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new ERICore(*this));
  }

  void update(const krims::ParameterMap& map) override {
    const std::string occ_coeff_key = Integral<stored_mtx_type>::update_key_coefficients;

    if (!map.exists(occ_coeff_key)) return;

    // Get coefficients as a shared pointer (having ownership)
    coefficients_occupied_ptr = static_cast<coefficients_ptr_type>(
          map.at_ptr<coefficients_type>(occ_coeff_key));

    // We will contract the coefficient row index over the number of
    // basis functions.
    if (coefficients_occupied_ptr->n_vectors() == 0) return;
    assert_size(coefficients_occupied_ptr->n_elem(),
                detail::Static14Data<stored_mtx_type>::nbas);
  }

  // TODO use static keys to generate values
  /** \brief Get the identifier of the integral */
  std::string id() const override {
    return std::string("static/14/ERI_") + (exchange ? "K" : "J");
  }

  /** \brief Get the friendly name of the integral */
  std::string name() const override {
    return std::string("Electron Repulsion Integrals, ") +
           (exchange ? "Exchange" : "Coulomb") + " operator";
  }

  ERICore(bool exchange, real_type k) : exchange(exchange), k(k) {}
};

}  // namespace static14
}  // namespace atomic
}  // namespace gint

#endif
