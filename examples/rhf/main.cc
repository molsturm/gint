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
#include <gint/IntegralUpdateKeys.hh>
#include <gint/Structure.hh>
#include <krims/GenMap.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/eigensystem.hh>

namespace rhf {
void restricted_hartree_fock(const krims::GenMap& intparams, const size_t n_alpha,
                             const size_t n_orbs = 4, const size_t max_iter = 100) {
  using gint::IntegralTypeKeys;
  using gint::IntegralUpdateKeys;
  using namespace lazyten;

  // Get the integrals for the appropriate integral parameters
  auto integrals = gint::IntegralLookup<SmallMatrix<double>>(intparams);
  auto S_bb      = integrals.lookup_integral(IntegralTypeKeys::overlap);
  auto T_bb      = integrals.lookup_integral(IntegralTypeKeys::kinetic);
  auto V_bb      = integrals.lookup_integral(IntegralTypeKeys::nuclear_attraction);
  auto J_bb      = integrals.lookup_integral(IntegralTypeKeys::coulomb);
  auto K_bb      = integrals.lookup_integral(IntegralTypeKeys::exchange);

  // Get an hcore guess
  const auto hcore_bb      = T_bb + V_bb;
  const auto eigensolution = eigensystem_hermitian(hcore_bb, S_bb, n_orbs);

  // Current occupied coefficients in convenient data structure
  const auto cocc = eigensolution.evectors().subview({0, n_alpha});

  // Initialise two-electron terms with guess coefficients
  J_bb.update({{IntegralUpdateKeys::coefficients_occupied, cocc}});
  K_bb.update({{IntegralUpdateKeys::coefficients_occupied, cocc}});

  double oldene = 0;
  std::cout << "Iter      etot      echange" << std::endl;
  for (size_t i = 0; i < max_iter; ++i) {
    // Obtain new eigenpairs ...
    const auto F_bb          = hcore_bb + (2 * J_bb - K_bb);
    const auto eigensolution = eigensystem_hermitian(F_bb, S_bb, n_orbs);

    // ... and new occupied coefficients
    const auto cocc = eigensolution.evectors().subview({0, n_alpha});

    // Compute HF energies
    // For example coulomb energy is tr(C^T J C)
    double ene_one_elec = trace(outer_prod_sum(cocc, hcore_bb * cocc));
    double ene_coulomb  = 2 * trace(outer_prod_sum(cocc, J_bb * cocc));
    double ene_exchge   = -trace(outer_prod_sum(cocc, K_bb * cocc));
    double energy       = 2 * (ene_one_elec + .5 * ene_coulomb + .5 * ene_exchge);

    // Display current iteration
    double energy_change = energy - oldene;
    std::cout << i << " " << energy << " " << energy_change << std::endl;
    oldene = energy;

    // Check for convergence
    if (fabs(energy_change) < 1e-6) break;

    // Update the two-electron integrals, before coefficients
    // go out of scope
    J_bb.update({{IntegralUpdateKeys::coefficients_occupied, cocc}});
    K_bb.update({{IntegralUpdateKeys::coefficients_occupied, cocc}});
  }

  std::cout << "Doubly occupied orbitals: " << std::endl;
  for (size_t i = 0; i < n_alpha; ++i) {
    std::cout << "  " << eigensolution.evalues()[i] << std::endl;
  }
}
}  // namespace rhf

int main() {
  gint::Structure be{
        {"Be", {{0., 0., 0.}}},
  };
  const size_t n_alpha = 2;

  krims::GenMap params_gauss{
        {"basis_type", "gaussian/libint"},
        {"basis_set_name", "pc-2"},
        {"structure", be},
  };

  krims::GenMap params_sturm{
        {"basis_type", "sturmian/atomic/cs_dummy"},
        {"k_exponent", 1.988},
        {"n_max", 5},
        {"l_max", 1},
        {"m_max", 1},
        {"structure", be},
  };

  rhf::restricted_hartree_fock(params_gauss, n_alpha);
  return 0;
}
