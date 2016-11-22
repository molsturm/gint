#pragma once
#ifdef GINT_STATIC_INTEGRALS
#include "Static14Data.hh"
#include "gint/Integral.hh"
#include "gint/IntegralCoreBase.hh"
#include <krims/ParameterMap.hh>

namespace gint {
namespace atomic {
namespace static14 {

template <typename StoredMatrix>
class OverlapIntegralCore;
template <typename StoredMatrix>
class NuclearAttractionIntegralCore;
template <typename StoredMatrix>
class KineticIntegralCore;
template <typename StoredMatrix, bool Exchange>
class ERICore;

/** This integral collection contains precomputed data for an atomic sturmian
 * basis with n_max = 3 and l_max=2, which yields 14 basis functions.
 *
 * The exponent k and the nuclear charge Z can still be chosen freely.
 */
template <typename StoredMatrix>
class IntegralCollection {
public:
  typedef StoredMatrix stored_matrix_type;
  typedef Integral<stored_matrix_type> integral_matrix_type;
  typedef typename stored_matrix_type::size_type size_type;
  typedef double real_type;

  static_assert(
        std::is_same<double, typename stored_matrix_type::scalar_type>::value,
        "The scalar type of StoredMatrix needs to be double for now.");

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
  integral_matrix_type operator()(const std::string& integral_id) const;
};

//
// Implementation
//
template <typename StoredMatrix>
const std::string IntegralCollection<StoredMatrix>::id = "atomic/static14";
template <typename StoredMatrix>
const std::string IntegralCollection<StoredMatrix>::name =
      "Fully precomputed integral data for 14 sturmians";

template <typename StoredMatrix>
IntegralCollection<StoredMatrix>::IntegralCollection(
      const krims::ParameterMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{parameters.at<double>("Z_charge")} {}

template <typename StoredMatrix>
Integral<StoredMatrix> IntegralCollection<StoredMatrix>::operator()(
      const std::string& integral_name) const {
  // TODO provide these keys as static constants somewhere
  if (integral_name == "nuclear_attraction")
    return make_integral<NuclearAttractionIntegralCore<stored_matrix_type>>(
          k_exponent, Z_charge);
  if (integral_name == "overlap")
    return make_integral<OverlapIntegralCore<stored_matrix_type>>();
  if (integral_name == "kinetic")
    return make_integral<KineticIntegralCore<stored_matrix_type>>(k_exponent);
  if (integral_name == "coulomb")
    return make_integral<ERICore<stored_matrix_type, false>>(k_exponent);
  if (integral_name == "exchange")
    return make_integral<ERICore<stored_matrix_type, true>>(k_exponent);

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<StoredMatrix>(nullptr);
}

//
// Integral cores
//

template <typename StoredMatrix>
class NuclearAttractionIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;

  const real_type k, Z;

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief return an element of the matrix \f$ {V_0}_{\mu',\mu} = -Zk/n
   * \delta_{\mu',\mu} \f$ */
  scalar_type operator()(size_type row, size_type col) const override {
    return -k * Z *
           detail::Static14Data<stored_matrix_type>::v0_bb_base(row, col);
  }

  bool has_transpose_operation_mode() const override {
    return detail::Static14Data<stored_matrix_type>::v0_bb_base
          .has_transpose_operation_mode();
  }

  void apply(const linalgwrap::MultiVector<
                   const linalgwrap::MutableMemoryVector_i<scalar_type>>& x,
             linalgwrap::MultiVector<
                   linalgwrap::MutableMemoryVector_i<scalar_type>>& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
             const scalar_type c_y =
                   linalgwrap::Constants<scalar_type>::zero) const override {
    detail::Static14Data<stored_matrix_type>::v0_bb_base.apply(
          x, y, mode, -k * Z * c_this, c_y);
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   */
  virtual void extract_block(
        stored_matrix_type& M, const size_type start_row,
        const size_type start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M =
              linalgwrap::Constants<scalar_type>::zero) const override {
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

template <typename StoredMatrix>
class OverlapIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief return an element of the matrix \f$ {S}_{\mu',\mu} */
  scalar_type operator()(size_type row, size_type col) const override {
    return detail::Static14Data<stored_matrix_type>::s_bb(row, col);
  }

  bool has_transpose_operation_mode() const override {
    return detail::Static14Data<stored_matrix_type>::s_bb
          .has_transpose_operation_mode();
  }

  void apply(const linalgwrap::MultiVector<
                   const linalgwrap::MutableMemoryVector_i<scalar_type>>& x,
             linalgwrap::MultiVector<
                   linalgwrap::MutableMemoryVector_i<scalar_type>>& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
             const scalar_type c_y =
                   linalgwrap::Constants<scalar_type>::zero) const override {
    detail::Static14Data<stored_matrix_type>::s_bb.apply(x, y, mode, c_this,
                                                         c_y);
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.  */
  virtual void extract_block(
        stored_matrix_type& M, const size_type start_row,
        const size_type start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M =
              linalgwrap::Constants<scalar_type>::zero) const override {
    detail::Static14Data<stored_matrix_type>::s_bb.extract_block(
          M, start_row, start_col, mode, c_this, c_M);
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

template <typename StoredMatrix>
class KineticIntegralCore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;

  const real_type k;

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief return an element of the matrix \f$ {T}_{\mu',\mu} */
  scalar_type operator()(size_type row, size_type col) const override {
    return k * k *
           detail::Static14Data<stored_matrix_type>::t_bb_base(row, col);
  }

  bool has_transpose_operation_mode() const override {
    return detail::Static14Data<stored_matrix_type>::t_bb_base
          .has_transpose_operation_mode();
  }

  void apply(const linalgwrap::MultiVector<
                   const linalgwrap::MutableMemoryVector_i<scalar_type>>& x,
             linalgwrap::MultiVector<
                   linalgwrap::MutableMemoryVector_i<scalar_type>>& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
             const scalar_type c_y =
                   linalgwrap::Constants<scalar_type>::zero) const override {
    detail::Static14Data<stored_matrix_type>::t_bb_base.apply(
          x, y, mode, k * k * c_this, c_y);
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.  */
  virtual void extract_block(
        stored_matrix_type& M, const size_type start_row,
        const size_type start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M =
              linalgwrap::Constants<scalar_type>::zero) const override {
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

template <typename StoredMatrix, bool Exchange>
class ERICore : public IntegralCoreBase<StoredMatrix> {
public:
  typedef IntegralCoreBase<StoredMatrix> base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef const linalgwrap::MultiVector<const vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  //! The exponent
  const real_type k;

  //! The occupied coefficients as a pointer
  coefficients_ptr_type coefficients_occupied_ptr;

  /** \brief Number of rows of the matrix */
  size_type n_rows() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override {
    return detail::Static14Data<stored_matrix_type>::nbas;
  }

  /** \brief return an element of the matrix \f$ {J}_{\mu',\mu} \f$ or \f$
   * {k}_{\mu',\mu} \f$. */
  scalar_type operator()(size_type a, size_type b) const override {
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
    const auto& i_bbbb = detail::Static14Data<stored_matrix_type>::i_bbbb_base;
    const size_type nbas = detail::Static14Data<stored_matrix_type>::nbas;

    assert_dbg(coefficients_occupied_ptr != nullptr,
               krims::ExcInvalidPointer());

    // Density matrix expression
    auto density_bb = outer_prod_sum(*coefficients_occupied_ptr,
                                     *coefficients_occupied_ptr);
    assert_dbg(density_bb.n_rows() == nbas, krims::ExcInternalError());
    assert_dbg(density_bb.n_cols() == nbas, krims::ExcInternalError());

    // Shell pair index for basis functions a and b:
    const size_type ab_pair = a * nbas + b;

    // Sum accumulator variable for this exchange or
    // coulomb  matrix element
    scalar_type mat_ab{0};

    // Double loop over basis functions c and d:
    for (size_type c = 0; c < nbas; ++c) {
      // Shell pair index for basis functions a and c:
      const size_type ac_pair = a * nbas + c;

      for (size_type d = 0; d < nbas; ++d) {
        // Shell pair index for basis functions c and d:
        // or b and d:
        const size_type cd_pair = c * nbas + d;
        const size_type db_pair = d * nbas + b;

        // Perform contraction:
        const scalar_type i_elem =
              Exchange ? i_bbbb(ac_pair, db_pair) : i_bbbb(ab_pair, cd_pair);
        mat_ab += density_bb(c, d) * k * i_elem;
      }  // d
    }    // c

    return mat_ab;
  }

  bool has_transpose_operation_mode() const override { return true; }

  void apply(const linalgwrap::MultiVector<
                   const linalgwrap::MutableMemoryVector_i<scalar_type>>& x,
             linalgwrap::MultiVector<
                   linalgwrap::MutableMemoryVector_i<scalar_type>>& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
             const scalar_type c_y =
                   linalgwrap::Constants<scalar_type>::zero) const override {
    using namespace linalgwrap;
    const size_type nbas = detail::Static14Data<stored_matrix_type>::nbas;

    assert_finite(c_this);
    assert_finite(c_y);
    assert_size(x.n_vectors(), y.n_vectors());
    assert_size(x.n_elem(), nbas);
    assert_size(y.n_elem(), nbas);
    assert_sufficiently_tested(mode != Transposed::ConjTrans);

    // Scale the current values of y or set them to zero
    // (if c_y == 0): We are now done with c_y and do not
    // need to worry about it any more in this function
    // TODO we are kind of calling an internal function here
    for (auto& vec : y) linalgwrap::detail::scale_or_set(vec, c_y);

    // if c_this == 0 we are done
    if (c_this == Constants<scalar_type>::zero) return;

    for (size_type veci = 0; veci < x.n_vectors(); ++veci) {
      const auto& vec = x[veci];
      auto& out = y[veci];
      for (size_type row = 0; row < nbas; ++row) {
        scalar_type sum = 0;
        for (size_type k = 0; k < nbas; ++k) {
          switch (mode) {
            case Transposed::None:
            case Transposed::Trans:
              // Same case since we are symmetric
              sum += (*this)(row, k) * vec(k);
              break;
            case Transposed::ConjTrans:
              // A variant of std::conj, which does not return a
              // complex data type if scalar is real only.
              // TODO we are kind of calling an internal function here
              linalgwrap::detail::ConjFctr mconj;
              sum += mconj((*this)(row, k)) * vec(k);
              break;
          }  // mode
        }    // k
        out(row) += c_this * sum;
      }  // row
    }    // veci
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.  */
  virtual void extract_block(
        stored_matrix_type& M, const size_type start_row,
        const size_type start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M =
              linalgwrap::Constants<scalar_type>::zero) const override {
    using namespace linalgwrap;
    const size_type nbas = detail::Static14Data<stored_matrix_type>::nbas;

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

    for (size_type row = 0; row < M.n_rows(); ++row) {
      for (size_type col = 0; col < M.n_cols(); ++col) {
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
            M(row, col) +=
                  c_this * mconj((*this)(start_col + col, start_row + row));
            break;
        }  // mode
      }    // col
    }      // row
  }

  std::unique_ptr<base_type> clone() const override {
    return std::unique_ptr<base_type>(new ERICore(*this));
  }

  void update(const krims::ParameterMap& map) override {
    const std::string occ_coeff_key =
          Integral<stored_matrix_type>::update_key_coefficients;

    if (!map.exists(occ_coeff_key)) return;

    // Get coefficients as a shared pointer (having ownership)
    coefficients_occupied_ptr = static_cast<coefficients_ptr_type>(
          map.at_ptr<coefficients_type>(occ_coeff_key));

    // We will contract the coefficient row index over the number of
    // basis functions.
    if (coefficients_occupied_ptr->n_vectors() == 0) return;
    assert_size(coefficients_occupied_ptr->n_elem(),
                detail::Static14Data<stored_matrix_type>::nbas);
  }

  // TODO use static keys to generate values
  /** \brief Get the identifier of the integral */
  std::string id() const override {
    return std::string("static/14/ERI_") + (Exchange ? "K" : "J");
  }

  /** \brief Get the friendly name of the integral */
  std::string name() const override {
    return std::string("Electron Repulsion Integrals, ") +
           (Exchange ? "Exchange" : "Coulomb") + " operator";
  }

  ERICore(real_type k) : k(k) {}
};

}  // namespace static14
}  // namespace atomic
}  // namespace gint

#endif
