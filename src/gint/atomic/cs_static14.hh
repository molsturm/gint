#pragma once
#ifdef GINT_STATIC_INTEGRALS

#include "Static14Data.hh"
#include "gint/Integral.hh"
#include "gint/IntegralCollectionBase.hh"
#include "gint/IntegralCoreBase.hh"
#include "gint/config.hh"

namespace gint {
namespace atomic {
namespace cs_static14 {

// In this namespace all things are real:
typedef real_type scalar_type;
typedef real_stored_mtx_type stored_mtx_type;
typedef real_multivector_type multivector_type;
typedef const_real_multivector_type const_multivector_type;

class OverlapIntegralCore;
class NuclearAttractionIntegralCore;
class KineticIntegralCore;
class ERICore;

void apply_stored_matrix(const real_stored_mtx_type& A, const_multivector_type& x,
                         multivector_type& y,
                         const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                         const scalar_type c_A = 1, const scalar_type c_y = 0);

/** This integral collection contains precomputed data for an atomic sturmian
 * basis with n_max = 3 and l_max=2, which yields 14 basis functions.
 *
 * The exponent k and the nuclear charge Z can still be chosen freely.
 */
class IntegralCollection final
      : public IntegralCollectionBase<OrbitalType::COMPLEX_ATOMIC> {
 public:
  typedef IntegralCollectionBase<OrbitalType::COMPLEX_ATOMIC> base_type;

  const static std::string id;
  const real_type k_exponent, Z_charge;

  /** Construct collection object from a set of parameters
   *
   * The following parameters are read:
   *   - k_exponent (double): The exponent of all Coulomb sturmians
   *   - Z_charge (double): The nuclear change of the system
   */
  IntegralCollection(const krims::GenMap& parameters);

  /** Lookup an integral by its type */
  integral_matrix_type lookup_integral(IntegralType type) const override;

  /** Obtain the id string of the collection / basis type */
  const std::string& basis_id() const override { return id; }

  /** Obtain the friendly name of the collection / basis type */
  std::string basis_name() const override {
    return "Fully precomputed integral data for 14 sturmians";
  }

  /** Create an integral collection for a particular basis set defined by parameters */
  static std::unique_ptr<base_type> create(const krims::GenMap& parameters) {
    return krims::make_unique<IntegralCollection>(parameters);
  }
};

//
// Integral cores
//
class NuclearAttractionIntegralCore final
      : public IntegralCoreBase<real_stored_mtx_type> {
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
  void extract_block(
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

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override {
    return {IntegralCollection::id, IntegralType::nuclear_attraction};
  }

  NuclearAttractionIntegralCore(real_type k, real_type Z) : k(k), Z(Z) {}
};

class OverlapIntegralCore final : public IntegralCoreBase<real_stored_mtx_type> {
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
  void extract_block(
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

  IntegralIdentifier id() const override {
    return {IntegralCollection::id, IntegralType::overlap};
  }
};

class KineticIntegralCore final : public IntegralCoreBase<real_stored_mtx_type> {
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
  void extract_block(
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

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override {
    return IntegralIdentifier{IntegralCollection::id, IntegralType::kinetic};
  }

  KineticIntegralCore(real_type k) : k(k) {}
};

class ERICore final : public IntegralCoreBase<real_stored_mtx_type> {
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
  scalar_type operator()(size_t a, size_t b) const override;

  bool has_transpose_operation_mode() const override { return true; }

  void apply(const_multivector_type& x, multivector_type& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_A = 1, const scalar_type c_y = 0) const override;

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.  */
  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
                     const scalar_type c_this = 1,
                     const scalar_type c_M = 0) const override;

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new ERICore(*this));
  }

  void update(const krims::GenMap& map) override;

  /** \brief Get the identifier of the integral */
  IntegralIdentifier id() const override {
    return {IntegralCollection::id,
            (exchange ? IntegralType::exchange : IntegralType::coulomb)};
  }

  ERICore(bool exchange, real_type k) : exchange(exchange), k(k) {}
};

}  // namespace cs_static14
}  // namespace atomic
}  // namespace gint

#endif
