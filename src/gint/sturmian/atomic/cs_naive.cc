#include "cs_naive.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_naive {

const std::string IntegralCollection::id = "atomic/cs_naive";

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : m_system{}, m_integral_calculator{} {
  const int n_max = parameters.at<int>("n_max");
  const int l_max = parameters.at<int>("l_max");
  const int m_max = parameters.at<int>("m_max");

  assert_throw(0 < n_max && n_max <= 6,
               ExcInvalidIntegralParameters(
                     "Maximum principle quantum number (" + std::to_string(n_max) +
                     ") needs to be in the range [1,6] for cs_naive, since higher "
                     "values are not yet implemented."));
  assert_throw(0 <= l_max && l_max <= 4,
               ExcInvalidIntegralParameters(
                     "Maximum angular momentum quantum number (" + std::to_string(l_max) +
                     ") needs to be in the range [0,4] for cs_naive, since higher values "
                     "are not yet implemented."));
  assert_throw(
        0 <= m_max && m_max <= 4,
        ExcInvalidIntegralParameters("Maximum magnetic momentum quantum number (" +
                                     std::to_string(m_max) +
                                     ") needs to be in the range [0,4] for cs_naive, "
                                     "since higher values are not yet implemented."));

  m_system.Z = parameters.at<scalar_type>("Z_charge");
  m_system.k = parameters.at<scalar_type>("k_exponent");
  m_system.basis = NlmBasis(n_max, l_max, m_max);
}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  typedef ERICore eri_type;
  const std::string& id = IntegralCollection::id;

  switch (type) {
    case IntegralType::nuclear_attraction:
      return make_integral<NuclearAttractionIntegralCore>(m_system, id);
    case IntegralType::overlap:
      return make_integral<OverlapIntegralCore>(m_system, id);
    case IntegralType::kinetic:
      return make_integral<KineticIntegralCore>(m_system, id);

    case IntegralType::coulomb: /* nobreak */
    case IntegralType::exchange:
      return make_integral<eri_type>(m_integral_calculator, m_system,
                                     IntegralIdentifier{id, type});

    default:
      assert_dbg(false, krims::ExcNotImplemented());
      return Integral<stored_matrix_type>(nullptr);
  }
}

//
// ERI Integral core
//

void ERICore::apply(const const_multivector_type& x, multivector_type& y,
                    const linalgwrap::Transposed mode, const scalar_type c_A,
                    const scalar_type c_y) const {
  // TODO Conceptionally this is code duplication with the apply function
  // in gint/IntegralCoreBase.cc, but just performed at elevated long double
  // precision. Be careful about this.
  //
  // If the elevated precision is not need, one might be able to remove this,
  // else one could merge it into the above place.

  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != linalgwrap::Transposed::ConjTrans);
  // All modes are same case since we are symmetric and real, so no
  // switching over mode.

  const size_t norb = x.n_rows(), n_vectors = x.n_cols();

  for (size_t i = 0; i < y.n_rows(); i++)
    for (size_t j = 0; j < y.n_cols(); j++) y(i, j) = (c_y != 0 ? c_y * y(i, j) : 0);

  // Work internally in highest available precision
  vector<working_scalar_type> ysum(n_vectors);
  for (size_t a = 0; a < norb; a++) {
    std::fill(std::begin(ysum), std::end(ysum), 0);
    for (size_t b = 0; b < norb; b++) {
      for (size_t q = 0; q < n_vectors; q++) {
        ysum[q] += c_A * compute_jk_element(a, b) * x(b, q);
      }
    }
    for (size_t q = 0; q < n_vectors; q++) y(a, q) += ysum[q];
  }
}

}  // namespace cs_naive
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
