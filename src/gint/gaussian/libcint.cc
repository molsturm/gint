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

#include "libcint.hh"

#ifdef GINT_HAVE_LIBCINT
#include "BasisSet.hh"
#include "IntegralLookupKeys.hh"
#include "Shell.hh"
#include "read_basisset.hh"
#include <krims/DataFiles/FindDataFile.hh>
#include <krims/FileSystem.hh>

// The interface of libcint is plain C, so the next section should be compiled with extern
// "C" linkage
extern "C" {
#include <cint.h>

//
// Unfortunately the cint.h header does not define all functions needed,
// so we need to add a few manually here.
//

/** The function which computes <i|j> inside libcint */
FINT cint1e_ovlp_sph(double* opij, const FINT* shls, const FINT* atm, const FINT natm,
                     const FINT* bas, const FINT nbas, const double* env);

/** The function which computes <i|V_nuc|j> inside libcint */
FINT cint1e_nuc_sph(double* opij, const FINT* shls, const FINT* atm, const FINT natm,
                    const FINT* bas, const FINT nbas, const double* env);

/** The function which computes <i| \delta j> inside libcint */
FINT cint1e_kin_sph(double* opij, const FINT* shls, const FINT* atm, const FINT natm,
                    const FINT* bas, const FINT nbas, const double* env);
}  // extern "C"

namespace gint {
namespace gaussian {
namespace libcint {

static_assert(std::is_same<FINT, int_type>::value,
              "The type used as int_type and the actual FINT integer type used by "
              "libcint need to agree.");

/** When initialised with a libcint::system it computes and makes available
 *  some data, which is needed to work with libcint shell pairs
 *
 *  TODO This exists in a similar way in libint. Perhaps one could create a
 *       generalisation which applies to both.
 */
struct BasisShellData {
  /** Return the index of the first basis function of a particular shell */
  size_t first_bfct(size_t shell) const { return m_bfct_range_idcs[shell]; }

  /** Return the number of basis functions of a particular shell */
  size_t n_bfct(size_t shell) const {
    return m_bfct_range_idcs[1 + shell] - m_bfct_range_idcs[shell];
  }

  /** Return the range of basis function indices covered by a particular shell */
  krims::Range<size_t> range_bfct(size_t shell) const {
    return {m_bfct_range_idcs[shell], m_bfct_range_idcs[1 + shell]};
  }

  /** Return the largest number of basis functions per shell */
  size_t max_n_bfct() const { return m_max_n_bfct; }

  //  /** Return the shell information for a particular shell */
  //  LibintShell shell_info(size_t shell) const {
  //    return LibintShell{shell, m_system_ptr->basis()[shell].size(), shell2bf[shell]};
  //  }

  /** Return the number of shells */
  size_t n_shells() const { return m_system_ptr->n_shells(); }

  /** Return the number of basis functions */
  size_t n_bas() const { return m_system_ptr->n_bas(); }

  /** Return the shell range which covers the requested index range
   * in the basis functions.
   *
   * In other words if only some basis function indices are required,
   * this function returns the range of shells which need to be computed
   * in order to compute the values for all basis functions.
   * Note that it is very likely that extra values at the beginning and
   * the end will be computed as well.
   */
  krims::Range<size_t> containing_shell_range(krims::Range<size_t> bfct_indices) const;

  explicit BasisShellData(const System& system);

 private:
  // Pointer to the molecular system info we use:
  krims::SubscriptionPointer<const System> m_system_ptr;

  /** The map shell index to index of first basis function of the shell.
   *  This map is one entry longer than the number of basis functions.
   *
   *  E.g. if we have 5 basis functions:
   *
   *       m_bfct_range_idcs[0] = index of the first basis function in shell 0
   *       m_bfct_range_idcs[1] = index of the first basis function in shell 1
   *       m_bfct_range_idcs[2] = index of the first basis function in shell 2
   *       m_bfct_range_idcs[3] = index of the first basis function in shell 3
   *       m_bfct_range_idcs[4] = index of the first basis function in shell 4
   *       m_bfct_range_idcs[5] = Total number of basis functions.
   */
  std::vector<size_t> m_bfct_range_idcs;

  // Maximal number of basis functions a shell may have
  size_t m_max_n_bfct;
};

BasisShellData::BasisShellData(const System& system)
      : m_system_ptr("BasisShellData", system),
        m_bfct_range_idcs(system.n_shells() + 1),
        m_max_n_bfct(0) {
  // Only this function is special to libcint!
  // If we abstract this we have the BasisShellData object be more general.
  size_t accu = 0;
  for (size_t i = 0; i <= system.n_shells(); ++i) {
    m_bfct_range_idcs[i] = accu;
    int_type n_bfct =
          CINTcgtos_spheric(static_cast<int_type>(i), system.shell_data.data());
    m_max_n_bfct = std::max(m_max_n_bfct, static_cast<size_t>(n_bfct));
    accu += static_cast<size_t>(n_bfct);
  }
  assert_internal(m_bfct_range_idcs[system.n_shells()] == system.n_bas());
}

krims::Range<size_t> BasisShellData::containing_shell_range(
      krims::Range<size_t> bfct_indices) const {
  // We assume here that shell2bf is sorted:
  assert_internal(std::is_sorted(m_bfct_range_idcs.begin(), m_bfct_range_idcs.end()));
  assert_internal(m_bfct_range_idcs.size() > 0);
  assert_internal(bfct_indices.upper_bound() <= m_system_ptr->n_bas());

  // Shortcut for the full range:
  if (bfct_indices.lower_bound() == 0 &&
      bfct_indices.upper_bound() == static_cast<size_t>(m_system_ptr->n_bas())) {
    return {0, n_shells()};
  }

  // Get the lower bound of the element one greater than the basis function index.
  // The lower bound is the first element which is not smaller than index+1,
  // so the one before that is the first shell we care about.
  const auto lower =
        std::lower_bound(std::begin(m_bfct_range_idcs), std::end(m_bfct_range_idcs),
                         bfct_indices.lower_bound() + 1);
  const size_t begin_shell =
        static_cast<size_t>(lower - std::begin(m_bfct_range_idcs)) - 1;

  // This may not happen since bfct_indices[n_shells()] == n_bas()
  assert_internal(lower != std::end(m_bfct_range_idcs));

  // Get the lower bound on the upper bound of the basis function indices.
  // This is the first shell we do not care about any more
  const auto upper =
        std::lower_bound(lower, std::end(m_bfct_range_idcs), bfct_indices.upper_bound());
  const size_t end_shell = static_cast<size_t>(upper - std::begin(m_bfct_range_idcs));

  // This may not happen since bfct_indices[n_shells()] == n_bas()
  assert_internal(upper != std::end(m_bfct_range_idcs));

  return {begin_shell, end_shell};
}

//
// System
//
System::System(krims::SubscriptionPointer<const Structure> structure_ptr, Basis basis)
      : shell_data(basis.n_shells() * BAS_SLOTS),
        atom_data(structure_ptr->n_atoms() * ATM_SLOTS),
        m_n_bas(0),
        m_n_shells(basis.n_shells()),
        m_n_atoms(structure_ptr->n_atoms()) {
  std::cerr << "libcint support is experimental. Some cases (like cc-pvdz on water) are "
               "known to produce the wrong results."
            << std::endl;

  // Libcint reserves the first entries of the env array for internal usage:
  env_back = env_back + PTR_ENV_START;

  assert_throw(env_back_off() + structure_ptr->n_atoms() * 3 < env.size(),
               ExcInvalidIntegralParameters(
                     "Cannot use libcint with this many atoms. You supplied us with a "
                     "structure consisting of " +
                     std::to_string(structure_ptr->n_atoms()) + " atoms."));

  // Copy the atoms and coordinates.
  for (size_t i = 0; i < structure_ptr->n_atoms(); ++i) {
    auto& atom = (*structure_ptr)[i];
    assert_throw(!atom.has_deviating_charge(),
                 ExcInvalidIntegralParameters(
                       "None of the atoms may have a charge different from the atomic "
                       "number for libcint. The problematic atom has atomic number " +
                       std::to_string(atom.atomic_number) + " but charge " +
                       std::to_string(atom.nuclear_charge)));

    atom_data[CHARGE_OF + ATM_SLOTS * i] = static_cast<int_type>(atom.atomic_number);

    // Copy atom coordinates and note where they are.
    atom_data[PTR_COORD + ATM_SLOTS * i] = static_cast<int_type>(env_back_off());
    env_back = std::copy(atom.coords.begin(), atom.coords.end(), env_back);

    // Force the use of point charges by default
    atom_data[NUC_MOD_OF + ATM_SLOTS * i] = POINT_NUC;
  }
  assert_internal(env_back_off() < env.size());

  // TODO libcint requires that all basis functions are corresponding to atoms
  //      therefore we need to make sure that no two atoms have the same coordinates
  for (size_t ai = 0; ai < structure_ptr->n_atoms(); ++ai) {
    for (size_t aj = ai + 1; aj < structure_ptr->n_atoms(); ++aj) {
      assert_throw((*structure_ptr)[ai].coords != (*structure_ptr)[aj].coords,
                   ExcInvalidIntegralParameters("Cannot use libcint if two atoms in the "
                                                "molecular structure have exactly the "
                                                "same set of coordinates (i.e. are "
                                                "located at exactly the same place."));
    }
  }  // atoms

  // Copy the basis set.
  for (size_t i = 0; i < basis.n_shells(); ++i) {
    auto& shell = basis[i];
    assert_throw(env_back_off() + shell.n_primitives() * 2 < env.size(),
                 ExcInvalidIntegralParameters("Cannot use libcint with this many basis "
                                              "functions since space in the env array is "
                                              "not sufficient."));

    assert_throw(shell.pure || shell.l < 2,
                 ExcInvalidIntegralParameters(
                       "libcint can only handle pure Gaussian shells at the moment"));

    size_t atom = 0;
    for (; atom < structure_ptr->n_atoms(); ++atom) {
      if ((*structure_ptr)[atom].coords == shell.origin) break;
    }
    assert_throw(atom < structure_ptr->n_atoms(),
                 ExcInvalidIntegralParameters(
                       "libcint can only be used if each shell corresponds to one atom. "
                       "Ghost shells which correspond to no real atom cannot be used."));

    assert_throw(shell.l < ANG_MAX,
                 ExcInvalidIntegralParameters(
                       "libcint only supports angular momentum up to " +
                       std::to_string(ANG_MAX) + ", but an angular momentum of " +
                       std::to_string(shell.l) + " was encountered."));
    shell_data[ANG_OF + BAS_SLOTS * i] = shell.l;

    shell_data[ATOM_OF + BAS_SLOTS * i]  = static_cast<int_type>(atom);
    shell_data[NPRIM_OF + BAS_SLOTS * i] = static_cast<int_type>(shell.n_primitives());

    // TODO
    // Here is some potential for optimisation.
    // Some basis sets like the Dunning ones use the same exponents
    // for multiple basis functions. This way we can make this
    // more clear by specifying nctr_of == 2 whilst still having the
    // same number primitives (saves time on evaluation)
    //
    // Similarly if we have the same *type* of at atom at multiple places in the
    // molecule. If we point to the same data section, we might enable some optimisation
    // tricks inside libcint and furthermore can save some storage.
    shell_data[NCTR_OF + BAS_SLOTS * i] = 1;

    assert_internal(shell.exponents.size() == shell.coefficients.size());
    shell_data[PTR_EXP + BAS_SLOTS * i] = static_cast<int_type>(env_back_off());
    env_back = std::copy(shell.exponents.begin(), shell.exponents.end(), env_back);

    // Libcint wants the coefficients to contain the normalisation coefficients of the
    // primitives, so we add this during the copying via a std::transform
    shell_data[PTR_COEFF + BAS_SLOTS * i] = static_cast<int_type>(env_back_off());
    env_back = std::transform(shell.coefficients.begin(), shell.coefficients.end(),
                              shell.exponents.begin(), env_back,
                              [&shell](double coeff, double exp) {
                                return coeff * CINTgto_norm(shell.l, exp);
                              });
  }  // shell_data

  // Cache the number of basis functions
  m_n_bas = static_cast<size_t>(
        CINTtot_cgto_spheric(shell_data.data(), static_cast<int_type>(basis.n_shells())));
}

//
// ERITensor
//
void ERITensor::compute_kernel(const std::array<krims::Range<size_t>, 4>& block,
                               kernel_type kernel) const {
  const System& system = *m_system_ptr;

  //
  // TODO Code duplication with the ERICores
  // TODO Parallelise here!
  //

  // Initialise the optimiser
  CINTOpt* opt = nullptr;
  cint2e_sph_optimizer(&opt, system.atom_data.data(),
                       static_cast<int_type>(system.n_atoms()), system.shell_data.data(),
                       static_cast<int_type>(system.n_shells()), system.env.data());

  const BasisShellData data{system};  // Cache some shell data
  const krims::Range<size_t> a_shells = data.containing_shell_range(block[0]);
  const krims::Range<size_t> b_shells = data.containing_shell_range(block[1]);
  const krims::Range<size_t> c_shells = data.containing_shell_range(block[2]);
  const krims::Range<size_t> d_shells = data.containing_shell_range(block[3]);

  // Allocate a temporary buffer which can hold the maximum number of elements
  //
  // Note: libcint uses the column-major, i.e. Fortran ordering, convention.
  // Therefore we need to be a little careful when using the values in
  // the values array returned.
  std::vector<double> values(data.max_n_bfct() * data.max_n_bfct() * data.max_n_bfct() *
                             data.max_n_bfct());

  // And another one which will be filled with a copy of the above values
  // but in row-major format.
  std::vector<double> values_rowmajor(values.size());

  for (size_t sa : a_shells) {
    for (size_t sb : b_shells) {
      for (size_t sc : c_shells) {
        for (size_t sd : d_shells) {
          // Call the function to get the integral values
          int_type shls[4] = {static_cast<int_type>(sa), static_cast<int_type>(sb),
                              static_cast<int_type>(sc), static_cast<int_type>(sd)};

          int_type all_zero = cint2e_sph(
                values.data(), shls, system.atom_data.data(),
                static_cast<int_type>(system.n_atoms()), system.shell_data.data(),
                static_cast<int_type>(system.n_shells()), system.env.data(), opt);

          // If all_zero is zero than we can skip calling the kernel.
          if (all_zero == 0) continue;

          // Copy them all into row-major:
          for (size_t a = 0; a < data.n_bfct(sa); ++a) {
            for (size_t b = 0; b < data.n_bfct(sb); ++b) {
              for (size_t c = 0; c < data.n_bfct(sc); ++c) {
                for (size_t d = 0; d < data.n_bfct(sd); ++d) {
                  const size_t i_ab = a * data.n_bfct(sb) + b;
                  const size_t i_abcd =
                        (i_ab * data.n_bfct(sc) + c) * data.n_bfct(sd) + d;

                  const size_t c_cd = c + data.n_bfct(sc) * d;
                  const size_t c_abcd =
                        a + data.n_bfct(sa) * (b + data.n_bfct(sb) * c_cd);

                  assert_internal(i_abcd < values_rowmajor.size());
                  assert_internal(c_abcd < values.size());
                  values_rowmajor[i_abcd] = values[c_abcd];
                }  // d
              }    // c
            }      // b
          }        // a

          kernel({{data.range_bfct(sa), data.range_bfct(sb), data.range_bfct(sc),
                   data.range_bfct(sd)}},
                 values_rowmajor.data());

        }  // sd
      }    // sc
    }      // sb
  }        // sa

  CINTdel_optimizer(&opt);
}

//
// LibintCollection
//
const std::string IntegralCollection::id = "gaussian/libcint";

Basis construct_basis(const krims::GenMap& parameters) {
  // TODO Generalise this: This function is pretty much the same for
  //      *all* gaussian basis set libraries.

  // The key: Either a file or the name of a basis set
  const std::string& key =
        parameters.at<const std::string>(IntegralLookupKeys::basis_set);

  const std::string errmsg =
        " Could not construct gaussian basis for the structure (provided via key \"" +
        IntegralLookupKeys::structure + "\") and the gaussian basis set \"" + key +
        "\" (provided via key \"" + IntegralLookupKeys::basis_set + "\")." +
        "\n\nDetails:\n--------\n";

  // If key is a valid file => Read it directly,
  // else interpret it as a basis set name and look it up in the library
  // of basis set files shipped with gint.
  const BasisSet basis_set = [&key, &errmsg] {
    try {
      return krims::path_exists(key) ? read_basisset(key) : lookup_basisset(key);
    } catch (const krims::ExcDatafileNotFound& e) {
      assert_throw(false, ExcInvalidIntegralParameters(errmsg + e.extra()));
    } catch (const ExcInvalidBasisSetFile& e) {
      assert_throw(false, ExcInvalidIntegralParameters(errmsg + e.extra()));
    }
    return BasisSet();
  }();

  auto structure_ptr = parameters.at_ptr<const Structure>(IntegralLookupKeys::structure);
  try {
    return Basis{*structure_ptr, basis_set};
  } catch (const ExcNoBasisForAtom& e) {
    assert_throw(false, ExcInvalidIntegralParameters(errmsg + e.extra()));
    return Basis{};
  }
}

System make_system(const krims::GenMap& parameters) {
  auto structure_ptr = parameters.at_ptr<const Structure>(IntegralLookupKeys::structure);

  // Check if an explicit basis is specified: Use it
  // else construct a basis from the parameters
  if (parameters.exists(IntegralLookupKeys::basis)) {
    return System(structure_ptr, parameters.at<const Basis>(IntegralLookupKeys::basis));
  } else {
    return System(structure_ptr, construct_basis(parameters));
  }
}

IntegralCollection::IntegralCollection(const krims::GenMap& parameters)
      : m_system{make_system(parameters)}, m_eri_tensor{m_system} {}

Integral<stored_matrix_type> IntegralCollection::lookup_integral(
      IntegralType type) const {
  switch (type) {
    case IntegralType::nuclear_attraction: /* nobreak */
    case IntegralType::overlap:            /* nobreak */
    case IntegralType::kinetic:
      return make_integral<OneElecIntegralCore>(type, m_system);

    case IntegralType::coulomb: /* nobreak */
    case IntegralType::exchange:
      return make_integral<ERICore>(type, m_system);
  }

  assert_dbg(false, krims::ExcNotImplemented());
  return Integral<stored_matrix_type>(nullptr);
}

//
// LibCintIntegralCoreBase
//
scalar_type LibCintIntegralCoreBase::operator()(size_t row, size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());
  const BasisShellData data{*m_system_ptr};

  scalar_type ret = 0;
  auto kernel = [&row, &col, &ret, &data](size_t s, size_t t, const scalar_type* values) {
    assert_internal(data.first_bfct(s) <= row);
    assert_internal(data.first_bfct(t) <= col);

    if (values == nullptr) return;  // All values are zero => ret will be zero, too

    // Find the one value we need and extract it.
    const size_t fs = row - data.first_bfct(s);
    const size_t ft = col - data.first_bfct(t);
    ret             = values[fs * data.n_bfct(t) + ft];
  };

  // Compute the one shell pair of the kernel we truly need:
  compute({row, row + 1}, {col, col + 1}, std::move(kernel));
  return ret;
}

void LibCintIntegralCoreBase::extract_block(stored_matrix_type& M, const size_t start_row,
                                            const size_t start_col,
                                            const lazyten::Transposed mode,
                                            const scalar_type c_A,
                                            const scalar_type c_M) const {
  assert_finite(c_A);
  assert_finite(c_M);
  assert_greater_equal(start_row + M.n_rows(), n_rows());
  assert_greater_equal(start_col + M.n_cols(), n_cols());
  assert_sufficiently_tested(mode != lazyten::Transposed::ConjTrans);

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  // Set elements of M to zero (if c_M == 0)
  // or scale them according to c_M.
  // This deals entirely with the coefficient c_M
  if (c_M == 0) {
    M.set_zero();
  } else {
    M *= c_M;
  }
  if (c_A == 0) return;  // If c_A == 0 we are done

  // The range of basis function indices we are interested in:
  const auto row_range = krims::range(start_row, start_row + M.n_rows());
  const auto col_range = krims::range(start_col, start_col + M.n_cols());

  const BasisShellData data{*m_system_ptr};
  auto kernel = [&row_range, &col_range, &M, &c_A, &data](size_t s, size_t t,
                                                          const scalar_type* values) {
    if (values == nullptr) return;  // The whole shell-pair was screened away

    // The first basis function indices we care about in this shell pair block:
    size_t i_beg = std::max(data.range_bfct(s).lower_bound(), row_range.lower_bound());
    size_t j_beg = std::max(data.range_bfct(t).lower_bound(), col_range.lower_bound());

    // The last basis function indices we a care about in this shell pair block
    const size_t i_end =
          std::min(data.range_bfct(s).upper_bound(), row_range.upper_bound());
    const size_t j_end =
          std::min(data.range_bfct(t).upper_bound(), col_range.upper_bound());

    // libcint evaluates batches of shells at once. Since a shell consists of
    // all harmonics for a radial parts we get a bunch of values.
    // The integrals are packed into shellsets in row-major (C-like) form,
    // so this iterates over the results and copy to the matrix M
    for (size_t i = i_beg; i < i_end; ++i) {
      for (size_t j = j_beg; j < j_end; ++j) {
        // Compute offset into values array:
        const size_t fs  = i - data.first_bfct(s);
        const size_t ft  = j - data.first_bfct(t);
        const size_t fst = fs * data.n_bfct(t) + ft;

        // Shift i and j for placement into M:
        const size_t i_shifted = i - row_range.lower_bound();
        const size_t j_shifted = j - col_range.lower_bound();
        assert_internal(fst < data.n_bfct(s) * data.n_bfct(t));

        M(i_shifted, j_shifted) += c_A * values[fst];
      }  // j
    }    // i
  };
  compute(row_range, col_range, std::move(kernel));
}

void LibCintIntegralCoreBase::apply(const const_multivector_type& x, multivector_type& y,
                                    const lazyten::Transposed mode, const scalar_type c_A,
                                    const scalar_type c_y) const {
  assert_finite(c_A);
  assert_finite(c_y);
  assert_size(x.n_cols(), y.n_cols());
  assert_size(x.n_rows(), n_cols());
  assert_size(y.n_rows(), n_rows());
  assert_sufficiently_tested(mode != lazyten::Transposed::ConjTrans);
  // All modes are same since we are symmetric and real, so no
  // switching over mode.

  // Scale the current values of out or set them to zero
  // (if c_y == 0): We are now done with c_y and do not
  // need to worry about it any more in this function
  if (c_y == 0) {
    y.set_zero();
  } else {
    y *= c_y;
  }
  if (c_A == 0) return;  // If c_A == 0 we are done

  const BasisShellData data{*m_system_ptr};
  auto kernel = [&x, &y, &c_A, &data](size_t s, size_t t, const scalar_type* values) {
    if (values == nullptr) return;  // The whole shell-pair was screened away

    // libcint evaluates batches of shells at once. Since a shell consists of
    // all harmonics for a radial parts we get a bunch of values.
    for (size_t fs = 0; fs < data.n_bfct(s); ++fs) {
      for (size_t ft = 0; ft < data.n_bfct(t); ++ft) {
        const size_t fst = fs * data.n_bfct(t) + ft;  // Index of function pair in result
        const size_t i   = data.first_bfct(s) + fs;   // Absolute index in first shell
        const size_t j   = data.first_bfct(t) + ft;   // Absolute index in second shell

        for (size_t vec = 0; vec < x.n_cols(); ++vec) {
          y(i, vec) += c_A * values[fst] * x(j, vec);
        }  // vec
      }    // j
    }      // i
  };
  compute(std::move(kernel));
}

//
// OneElecIntegralCore
//
OneElecIntegralCore::OneElecIntegralCore(IntegralType type, const System& system)
      : base_type(system), m_cint1e_kernel(nullptr), m_type(type) {
  switch (m_type) {
    case IntegralType::overlap:
      m_cint1e_kernel = &cint1e_ovlp_sph;
      break;
    case IntegralType::nuclear_attraction:
      m_cint1e_kernel = &cint1e_nuc_sph;
      break;
    case IntegralType::kinetic:
      m_cint1e_kernel = &cint1e_kin_sph;
      break;
    default:
      assert_internal(false);
  }
}

void OneElecIntegralCore::compute(const krims::Range<size_t>& rows,
                                  const krims::Range<size_t>& cols,
                                  typename base_type::kernel_type&& kernel) const {
  const auto& system = base_type::system();
  const BasisShellData data{system};
  const krims::Range<size_t> row_shells = data.containing_shell_range(rows);
  const krims::Range<size_t> col_shells = data.containing_shell_range(cols);

  //
  // TODO Parallelise here!
  //

  // Allocate a temporary buffer which can hold the maximum number
  std::vector<double> values(data.max_n_bfct() * data.max_n_bfct() * data.max_n_bfct() *
                             data.max_n_bfct());

  // Loop over shell sets {s,t}
  for (size_t sa : row_shells) {
    for (size_t sb : col_shells) {
      // Compute the integral for a shell pair and call the kernel function with the
      // location where the computed integrals are stored.
      //
      // Note: libcint uses the column-major, i.e. Fortran ordering, convention
      // Since the one-electron integrals are all symmetric, we do not have to
      // care about this here, though
      int_type shls[2]  = {static_cast<int_type>(sa), static_cast<int_type>(sb)};
      int_type all_zero = m_cint1e_kernel(
            values.data(), shls, system.atom_data.data(),
            static_cast<int_type>(system.n_atoms()), system.shell_data.data(),
            static_cast<int_type>(system.n_shells()), system.env.data());

      // If all_zero == 0, then the values are all zero, hence we pass a nullptr
      kernel(sa, sb, all_zero == 0 ? nullptr : values.data());
    }  // sb
  }    // sa
}

//
// ERI Integral core
//
void ERICore::compute(const krims::Range<size_t>& rows, const krims::Range<size_t>& cols,
                      typename base_type::kernel_type&& kernel) const {
  const System& system(base_type::system());

  // Initialise the optimiser
  CINTOpt* opt = nullptr;
  cint2e_sph_optimizer(&opt, system.atom_data.data(),
                       static_cast<int_type>(system.n_atoms()), system.shell_data.data(),
                       static_cast<int_type>(system.n_shells()), system.env.data());

  //
  // TODO Parallelise here!
  // TODO Code duplication with the ERITensor
  //

  const BasisShellData data{system};
  const krims::Range<size_t> row_shells = data.containing_shell_range(rows);
  const krims::Range<size_t> col_shells = data.containing_shell_range(cols);

  // Allocate a temporary buffer which can hold the maximum number
  std::vector<double> values(data.max_n_bfct() * data.max_n_bfct() * data.max_n_bfct() *
                             data.max_n_bfct());

  // Compute density matrix:
  const auto dens = compute_density_matrix();

  // Loop over shell sets {sa,sb,sc,sd}
  for (size_t sa : row_shells) {
    for (size_t sb : col_shells) {
      // Allocate a buffer to contain the values in the shell pair {sa,sb}
      // (stored in row-major)
      std::vector<scalar_type> buffer(data.n_bfct(sa) * data.n_bfct(sb), 0.);

      // Flag to indicate that all computed values of the inner double loop are zero
      bool all_shellpairs_zero = true;

      for (size_t sc = 0; sc < data.n_shells(); ++sc) {
        for (size_t sd = 0; sd < data.n_shells(); ++sd) {
          // Swap shell indices sa and sc if computing exchange
          const bool exchange = m_type == IntegralType::exchange;
          size_t sA           = exchange ? sc : sa;
          size_t sC           = exchange ? sa : sc;

          // Call the function to get the integral values
          //
          // Note: libcint uses the column-major, i.e. Fortran ordering, convention.
          // Therefore we need to be a little careful when using the values in
          // the values array returned.
          int_type shls[4] = {static_cast<int_type>(sA), static_cast<int_type>(sb),
                              static_cast<int_type>(sC), static_cast<int_type>(sd)};
          int_type all_zero = cint2e_sph(
                values.data(), shls, system.atom_data.data(),
                static_cast<int_type>(system.n_atoms()), system.shell_data.data(),
                static_cast<int_type>(system.n_shells()), system.env.data(), opt);

          // If all_zero is zero than we can skip calling the kernel.
          if (all_zero == 0) continue;

          all_shellpairs_zero = false;
          for (size_t a = 0; a < data.n_bfct(sa); ++a) {
            for (size_t b = 0; b < data.n_bfct(sb); ++b) {
              // row-major index into buffer
              const size_t i_ab = a * data.n_bfct(sb) + b;

              for (size_t c = 0; c < data.n_bfct(sc); ++c) {
                // Swap shell indices a and c for exchange
                const size_t A = exchange ? c : a;
                const size_t C = exchange ? a : c;

                // Compute partial offsets (in column major form,
                // since values is a column-major, Fortran-like array)
                const size_t c_bC  = b + data.n_bfct(sb) * C;
                const size_t c_AbC = A + data.n_bfct(sA) * c_bC;
                const size_t d_stride =
                      data.n_bfct(sA) * data.n_bfct(sb) * data.n_bfct(sC);

                for (size_t d = 0; d < data.n_bfct(sd); ++d) {
                  // Indices to access element dens(c,d) of the density matrix
                  const size_t idens_c = data.first_bfct(sc) + c;
                  const size_t idens_d = data.first_bfct(sd) + d;

                  // Compute index into integral values (column-major!)
                  const size_t c_AbCd = c_AbC + d_stride * d;

                  // Perform explicit bound checks:
                  assert_internal(i_ab < buffer.size());
                  assert_internal(c_AbCd < data.n_bfct(sa) * data.n_bfct(sb) *
                                                 data.n_bfct(sc) * data.n_bfct(sd));

                  buffer[i_ab] += dens(idens_c, idens_d) * values[c_AbCd];
                }  // d
              }    // c
            }      // b
          }        // a

        }  // sh_d
      }    // sh_c

      // If all integrals have returned as zero, then pass only a nullptr
      kernel(sa, sb, all_shellpairs_zero ? nullptr : buffer.data());
    }  // sh_b
  }    // sh_a

  CINTdel_optimizer(&opt);
}

}  // namespace libcint
}  // namespace gaussian
}  // namespace gint
#endif  // GINT_HAVE_LIBCINT
