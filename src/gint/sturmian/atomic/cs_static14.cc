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

#include "cs_static14.hh"

#ifdef GINT_HAVE_STATIC_INTEGRALS
#include "gint/IntegralLookupKeys.hh"
#include "gint/IntegralUpdateKeys.hh"
#include "gint/Structure.hh"
#include "gint/sturmian/atomic/NlmBasis.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace cs_static14 {

const std::string IntegralCollection::id = "sturmian/atomic/cs_static14";

double atom_charge(const krims::GenMap& parameters) {
  if (!parameters.exists(IntegralLookupKeys::structure)) {
    // TODO Is this a sensible way to deal with this?
    // No structure provided, i.e. no atom, i.e. charge is zero
    return 0;
  }
  const auto structure_ptr =
        parameters.at_ptr<const Structure>(IntegralLookupKeys::structure);

  if (structure_ptr->n_atoms() == 0) return 0;
  assert_throw(structure_ptr->n_atoms() == 1,
               ExcInvalidIntegralParameters("The structure provided to cs_static14 is "
                                            "not an atom, but a molecule consisting of " +
                                            std::to_string(structure_ptr->n_atoms()) +
                                            " atoms."));
  return (*structure_ptr)[0].nuclear_charge;
}

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : k_exponent{parameters.at<double>("k_exponent")},
        Z_charge{atom_charge(parameters)},
        m_eri_tensor(k_exponent) {

  if (parameters.exists("nlm_basis")) {
    const auto& nlmbasis = parameters.at<const NlmBasis>("nlm_basis");

    // Construct an NlmBasis in normal, dense form with nlm order
    const int n_max = nlmbasis.n_max();
    const int l_max = nlmbasis.l_max();
    const int m_max = nlmbasis.m_max();
    const NlmBasis nlmbasis_compare(n_max, l_max, m_max);

    // Assert its the same
    assert_throw(nlmbasis_compare == nlmbasis,
                 ExcInvalidIntegralParameters("cs_static14 only works with an nlm_basis "
                                              "which is in nlm order and has no gaps up "
                                              "to (n_max,l_max,m_max) == (3,2,2)"));

  } else {
    const int n_max = parameters.at<int>("n_max");
    const int l_max = parameters.at<int>("l_max");
    const int m_max = parameters.at<int>("m_max");

    assert_throw(
          n_max == 3 && l_max == 2 && m_max == 2,
          ExcInvalidIntegralParameters(
                "cs_static14 only works with (n_max,l_max,m_max) == (3,2,2) and not (" +
                std::to_string(n_max) + "," + std::to_string(l_max) + "," +
                std::to_string(m_max) + ")"));
  }
}

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
      assert_implemented(false);
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
      assert_internal(false);
  }
}

void apply_stored_matrix(const stored_matrix_type& A, const const_multivector_type& x,
                         multivector_type& y, const lazyten::Transposed mode,
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
void ERITensor::compute_kernel(const std::array<krims::Range<size_t>, 4>& block,
                               kernel_type kernel) const {
  const auto& i_bbbb = Static14Data::i_bbbb_base;
  const size_t nbas  = Static14Data::nbas;

  // Copy the tensor data into a buffer and call the kernel
  std::vector<scalar_type> buffer(block[0].length() * block[1].length() *
                                  block[2].length() * block[3].length());

  auto it = std::begin(buffer);
  for (auto a : block[0]) {
    for (auto b : block[1]) {
      for (auto c : block[2]) {
        for (auto d : block[3]) {
          const size_t ab = a * nbas + b;
          const size_t cd = c * nbas + d;

          *(it++) = k * i_bbbb(ab, cd);
        }  // d
      }    // c
    }      // b
  }        // a

  kernel(block, buffer.data());
}

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
  const size_t nbas  = Static14Data::nbas;
  // Density matrix expression
  const auto density_bb = compute_density_matrix();
  assert_internal(density_bb.n_rows() == nbas && density_bb.n_cols() == nbas);

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
