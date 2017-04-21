#ifdef GINT_HAVE_STATIC_INTEGRALS
#include "cs_static14.hh"
#include "gint/IntegralLookupKeys.hh"
#include "gint/IntegralUpdateKeys.hh"
#include "gint/Structure.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_static14 {

const std::string IntegralCollection::id = "atomic/cs_static14";

double atom_charge(const krims::GenMap& parameters) {
  const auto structure_ptr =
        parameters.at_ptr<const Structure>(IntegralLookupKeys::structure);
  assert_throw(structure_ptr->n_atoms(),
               ExcInvalidIntegralParameters("The structure provided to cs_static14 is "
                                            "not an atom, but a molecule consisting of " +
                                            std::to_string(structure_ptr->n_atoms()) +
                                            " atoms."));
  return (*structure_ptr)[0].nuclear_charge;
}

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{atom_charge(parameters)} {}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  switch (type) {
    case IntegralType::nuclear_attraction: /* nobreak */
    case IntegralType::kinetic:            /* nobreak */
    case IntegralType::overlap:
      return make_integral<OneElecIntegralCore>(type, Z_charge, k_exponent);

    case IntegralType::coulomb: /* nobreak */
    case IntegralType::exchange:
      return make_integral<ERICore>(type, k_exponent);

    default:
      assert_throw(false, krims::ExcNotImplemented());
      return Integral<stored_matrix_type>(nullptr);
  }
}

OneElecIntegralCore::OneElecIntegralCore(IntegralType type, const scalar_type Z,
                                         const scalar_type k)
      : m_fac(1),
        m_mat_ptr("OneElecIntegralCore"),
        m_inv_mat_ptr("OneElecIntegralCore"),
        m_type(type) {
  switch (m_type) {
    case IntegralType::nuclear_attraction:
      m_fac = -k * Z;
      m_mat_ptr.reset(Static14Data::v0_bb_base);
      break;
    case IntegralType::overlap:
      m_mat_ptr.reset(Static14Data::s_bb);
      m_inv_mat_ptr.reset(Static14Data::sinv_bb);
      break;
    case IntegralType::kinetic:
      m_fac = k * k;
      m_mat_ptr.reset(Static14Data::t_bb_base);
      break;
    default:
      assert_dbg(false, krims::ExcInternalError());
  }
}

void apply_stored_matrix(const stored_matrix_type& A, const const_multivector_type& x,
                         multivector_type& y, const linalgwrap::Transposed mode,
                         const scalar_type c_A, const scalar_type c_y) {
  // scale y by c_y or set to zero
  for (size_t i = 0; i < x.n_rows(); i++)
    for (size_t j = 0; j < x.n_cols(); j++) y(i, j) = (c_y == 0 ? 0 : c_y * y(i, j));

  // Everything in this module is real and symmetric, so we can ignore mode.
  (void)mode;
  for (size_t k = 0; k < x.n_cols(); k++)  // Iterate over vectors
    for (size_t i = 0; i < A.n_rows(); i++) {
      real_type row_sum = 0;
      for (size_t j = 0; j < A.n_cols(); j++) row_sum += c_A * A(i, j) * x(j, k);
      y(i, k) += row_sum;
    }
}

//
// IntegralCores
//

typename ERICore::scalar_type ERICore::operator()(size_t a, size_t b) const {
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
  const auto& i_bbbb = Static14Data::i_bbbb_base;
  const size_t nbas = Static14Data::nbas;
  // Density matrix expression
  const auto density_bb = compute_density_matrix();
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
      const scalar_type i_elem = type == IntegralType::exchange
                                       ? i_bbbb(ac_pair, db_pair)
                                       : i_bbbb(ab_pair, cd_pair);
      mat_ab += k * i_elem * density_bb(c, d);
    }  // d
  }    // c

  return mat_ab;
}

}  // namespace cs_static14
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint

#endif  // ifdef GINT_HAVE_STATIC_INTEGRALS
