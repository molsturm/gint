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

#pragma once
#include <krims/ExceptionSystem.hh>

namespace gint {

/** Quantum number utility namespace */
namespace qnum {
DefException1(ExcInvalidAmLetter, char, << "Invalid angular momentum letter: " << arg1
                                        << ".");

/** Converts an AM letter to the am quantum number.
 *
 * \throws ExcInvalidAmLetter in case the letter is not a valid angular momentum letter.
 */
int letter_to_am(char l);

/** Converts an integer to an angular momentum quantum letter in upper case.
 *
 * \throws ExcOutsideRange in case the value of l is outside range
 */
char am_to_letter(int l);

}  // namespace qnum
}  // namespace gint
