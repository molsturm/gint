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

#include <gint/IntegralLookup.hh>
#include <krims/GenMap.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/eigensystem.hh>

namespace laplace {
typedef double scalar_type;
typedef lazyten::SmallMatrix<scalar_type> matrix_type;
typedef lazyten::SmallVector<scalar_type> vector_type;

auto spectrum(const krims::GenMap& intparams, const size_t n_ep = 10,
              const krims::GenMap& solveparams = krims::GenMap{})
      -> lazyten::Eigensolution<double, vector_type> {
  using gint::IntegralTypeKeys;

  // Get the kinetic and overlap integrals for the appropriate integral parameters
  auto integrals = gint::IntegralLookup<matrix_type>(intparams);
  auto S_bb      = integrals.lookup_integral(IntegralTypeKeys::overlap);
  auto T_bb      = integrals.lookup_integral(IntegralTypeKeys::kinetic);

  // Solve the eigensystem and return the result
  return lazyten::eigensystem_hermitian(T_bb, S_bb, n_ep, solveparams);
}
}  // namespace laplace

template <typename Solution>
void print_spectrum(const Solution& solution) {
  std::cout << "Obtained Laplace eigenvalues: " << std::endl;
  for (const auto& val : solution.evalues()) {
    std::cout << "  " << val << std::endl;
  }
}

int main() {
  {  // Solve using a small Coulomb-Sturmian basis.
    krims::GenMap params{
          {"basis_type", "sturmian/atomic/cs_static14"},
          {"k_exponent", 1.5},
          {"n_max", 3},
          {"l_max", 2},
          {"m_max", 2},
    };
    auto solution = laplace::spectrum(params);
    print_spectrum(solution);
  }

  return 0;
}
