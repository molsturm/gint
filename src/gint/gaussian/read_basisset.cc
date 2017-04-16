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

#include "read_basisset.hh"
#include "BasisSet.hh"
#include "Shell.hh"
#include "gint/Element.hh"
#include <algorithm>
#include <functional>

namespace gint {
namespace gaussian {

namespace detail {

// TODO Make publicly accessible utility function!
/** Converts an AM letter to the am quantum number. Returns -1 on error */
int letter_to_am(char l) {
  switch (std::toupper(l)) {
    /* clang-format off */
    case 'S': return 0;
    case 'P': return 1;
    case 'D': return 2;
    case 'F': return 3;
    case 'G': return 4;
    case 'H': return 5;
    case 'I': return 6;
    case 'J': return 7;
    case 'K': return 8;
    case 'L': return 9;
    case 'M': return 10;
    case 'N': return 11;
    case 'O': return 12;
    /* clang-format on */
    default:
      return -1;
  }
}

/** Convert a Fortran float string (containing 'D's and 'd's to a corresponding string
 * which can be interpreted by C++ */
std::string floatstring_fortran_to_c(std::string in) {
  std::for_each(std::begin(in), std::end(in), [](char& c) {
    if (c == 'd') c = 'e';
    if (c == 'D') c = 'E';
  });
  return in;
}

std::string to_format_string(BasisSetFileFormat fmt) {
  switch (fmt) {
    case BasisSetFileFormat::Gaussian94:
      return "Gaussian94";
    default:
      return "<unknown>";
  }
  return "<unknown>";
}

}  // namespace detail

//
// Gaussian94
//
template <size_t N>
void read_gaussian94_shell_section(
      std::istream& f, std::vector<real_type>& exponents,
      std::array<std::reference_wrapper<std::vector<real_type>>, N> coeff_vecs) {
// This function deals both with the "normal" cases of a single angular momentum
// shell(N=1), as well as the SP case(N=2)  where two coefficients are supplied at once
#ifdef DEBUG
  for (std::vector<real_type>& cv : coeff_vecs) {
    assert_dbg(exponents.size() == cv.size(), krims::ExcInternalError());
  }
#endif  // DEBUG

  std::string line;
  size_t i = 0;
  for (; i < exponents.size() && getline(f, line); ++i) {
    if (line.empty() || line[0] == '!') continue;

    std::stringstream ss(detail::floatstring_fortran_to_c(line));
    ss >> exponents[i];
    for (std::vector<real_type>& cv : coeff_vecs) ss >> cv[i];
    assert_throw(
          ss, ExcInvalidBasisSetFile("Could not read exponent and/or coefficient for " +
                                     std::to_string(i + 1) + "-th primitive."));
  }

  assert_throw(i == exponents.size(),
               ExcInvalidBasisSetFile("Number of expected primitives (" +
                                      std::to_string(exponents.size()) +
                                      ") and number of obtained primitives (" +
                                      std::to_string(i) + ") differs."));
}

std::vector<Shell> read_gaussian94_element_section(std::istream& f) {
  std::vector<Shell> ret;

  // Should the use of Cartesian d functions be enforced
  // By default this is not done, i.e. spherical harmonics are used instead.
  const bool force_cartesian_d = false;

  std::string line;
  while (std::getline(f, line) && line != "****") {
    if (line.empty() || line[0] == '!') continue;  // Ignore comments, empty

    std::string amstr;
    size_t n_contr;
    std::stringstream(line) >> amstr >> n_contr;

    assert_throw(n_contr >= 1,
                 ExcInvalidBasisSetFile("Invalid number of contractions encountered."));

    if (amstr == "sp" || amstr == "SP") {
      // This is a special SP shell, which we will represent as a
      // separate s and a separate p shell

      // S and P shell with no exponents/coefficients, but centered at the origin
      Shell s{0, false, {}, {}, {{0, 0, 0}}};
      Shell p{1, false, {}, {}, {{0, 0, 0}}};

      s.exponents.resize(n_contr);
      s.coefficients.resize(n_contr);
      p.coefficients.resize(n_contr);

      try {
        read_gaussian94_shell_section<2>(
              f, s.exponents, {{std::ref(s.coefficients), std::ref(p.coefficients)}});
      } catch (ExcInvalidBasisSetFile& e) {
        e.prepend_extra("Error when parsing " + std::to_string(1 + ret.size()) +
                        "-th shell with angular momentum string \"" + amstr + "\": ");
        throw;
      }

      p.exponents = s.exponents;
      ret.push_back(std::move(s));
      ret.push_back(std::move(p));
    } else {
      Shell sh;
      sh.l = detail::letter_to_am(amstr[0]);
      assert_throw(
            amstr.size() == 1 && sh.l != -1,
            ExcInvalidBasisSetFile("Invalid angular momentum string \"" + amstr + "\"."));

      sh.pure = force_cartesian_d ? (sh.l > 2) : (sh.l > 1);
      sh.exponents.resize(n_contr);
      sh.coefficients.resize(n_contr);
      sh.origin = {{0, 0, 0}};

      try {
        read_gaussian94_shell_section<1>(f, sh.exponents, {{std::ref(sh.coefficients)}});
      } catch (ExcInvalidBasisSetFile& e) {
        e.prepend_extra("Error when parsing " + std::to_string(1 + ret.size()) +
                        "-th shell with angular momentum string \"" + amstr + "\".");
        throw;
      }
      ret.push_back(std::move(sh));
    }
  }

  return ret;
}

BasisSet read_gaussian94(std::istream& f) {
  /** The format is roughly speaking:
   *
   *   ****
   *   Element_symbol
   *   AM   n_contr
   *      exp1     coeff1
   *      exp2     coeff2
   *   AM2  n_contr
   *      exp1     coeff1
   *      exp2     coeff2
   *   ****
   *   Element_symbol
   *   ...
   *
   * where
   *    Element_symbol    a well-known element symbol
   *    AM                azimuthal quantum number as a symbol for angular momentum
   *    n_contr           number of contracted primitives
   *    exp1              exponents
   *    coeff1            corresponding contraction coefficients.
   *
   * All element blocks are separated by '****'. The numbers for exp and coeff
   * may be given in the Fortran float convention (using 'D's instead of 'E's).
   */

  assert_throw(f, krims::ExcIO());

  // Seek until first "****":
  std::string line;
  while (std::getline(f, line) && line != "****") {
  }
  assert_throw(f,
               ExcInvalidBasisSetFile("EOF while seeking for initial '****' sequence"));

  BasisSet b;
  while (std::getline(f, line)) {
    // Skip empty lines or comment lines:
    if (line.empty() || line[0] == '!') continue;

    // Determine atomic number from element symbol
    std::stringstream ss(line);
    std::string symbol;
    ss >> symbol;

    assert_throw(is_element_symbol(symbol),
                 ExcInvalidBasisSetFile("Unknown element symbol: \"" + symbol + "\"."));
    const unsigned int at_num = Element::by_symbol(symbol).atomic_number;

    // Parse all shells for this element
    try {
      b.atomic_number_to_shells[at_num] = read_gaussian94_element_section(f);
    } catch (ExcInvalidBasisSetFile& e) {
      e.prepend_extra("Error in section corresponding to element \"" + symbol + "\": ");
      throw;
    }
  }

  // TODO How about normalisation ??
  //      Check how libint does this!

  return b;
}

//
// General
//

BasisSet read_basisset(std::istream& in, BasisSetFileFormat fmt) {
  try {
    switch (fmt) {
      case BasisSetFileFormat::Gaussian94:
        return read_gaussian94(in);
        break;
      default:
        assert_throw(false, krims::ExcNotImplemented());
    }
  } catch (ExcInvalidBasisSetFile& e) {
    e.prepend_extra(
          "The input stream provided to read_basisset does not contain a valid " +
          detail::to_format_string(fmt) + " file: ");
    throw;
  }
}

}  // namespace gaussian
}  // namespace gint
