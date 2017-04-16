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

#include "IntegralCollectionBase.hh"

namespace gint {
namespace sturmian {
namespace atomic {
namespace nlm_order {

const std::string IntegralCollectionBaseKeys::Z_charge = "Z_charge";
const std::string IntegralCollectionBaseKeys::k_exponent = "k_exponent";
const std::string IntegralCollectionBaseKeys::nlm_basis = "nlm_basis";

IntegralCollectionBase::IntegralCollectionBase(const krims::GenMap& parameters) {
  const auto Z = parameters.at<scalar_type>(IntegralCollectionBaseKeys::Z_charge);
  const auto k = parameters.at<scalar_type>(IntegralCollectionBaseKeys::k_exponent);

  if (parameters.exists(IntegralCollectionBaseKeys::nlm_basis)) {
    const auto& nlmbasis =
          parameters.at<const NlmBasis>(IntegralCollectionBaseKeys::nlm_basis);

    m_system = SturmintSystem(Z, k, nlmbasis);
  } else {
    const int n_max = parameters.at<int>("n_max");
    const int l_max = parameters.at<int>("l_max");
    const int m_max = parameters.at<int>("m_max");

    assert_throw(0 < n_max, ExcInvalidIntegralParameters(
                                  "Principle quantum number (" + std::to_string(n_max) +
                                  ") needs to be greater than zero."));
    assert_throw(
          0 <= l_max && l_max < n_max,
          ExcInvalidIntegralParameters("Maximum angular momentum quantum number (" +
                                       std::to_string(l_max) +
                                       ") needs to be in the range [0,n_max] == [0," +
                                       std::to_string(n_max) + "]."));
    assert_throw(
          0 <= m_max && m_max <= l_max,
          ExcInvalidIntegralParameters("Maximum magnetic momentum quantum number (" +
                                       std::to_string(m_max) +
                                       ") needs to be in the range [0,l_max] == [0," +
                                       std::to_string(l_max) + "]."));

    m_system = SturmintSystem(Z, k, NlmBasis(n_max, l_max, m_max));
  }
}

}  // namespace nlm_order
}  // namespace atomic
}  // namespace sturmian
}  // namespace gint
