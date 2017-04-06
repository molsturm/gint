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

#include "Element.hh"
#include <algorithm>

namespace gint {

const std::array<Element, 118>& elements() {
  static std::array<Element, 118> elements{{
        {1, "H", "hydrogen"},       {2, "He", "helium"},
        {3, "Li", "lithium"},       {4, "Be", "beryllium"},
        {5, "B", "boron"},          {6, "C", "carbon"},
        {7, "N", "nitrogen"},       {8, "O", "oxygen"},
        {9, "F", "fluorine"},       {10, "Ne", "neon"},
        {11, "Na", "sodium"},       {12, "Mg", "magnesium"},
        {13, "Al", "aluminium"},    {14, "Si", "silicon"},
        {15, "P", "phosphorus"},    {16, "S", "sulphur"},
        {17, "Cl", "chlorine"},     {18, "Ar", "argon"},
        {19, "K", "potassium"},     {20, "Ca", "calcium"},
        {21, "Sc", "scandium"},     {22, "Ti", "titanium"},
        {23, "V", "vanadium"},      {24, "Cr", "chromium"},
        {25, "Mn", "manganese"},    {26, "Fe", "iron"},
        {27, "Co", "cobalt"},       {28, "Ni", "nickel"},
        {29, "Cu", "copper"},       {30, "Zn", "zinc"},
        {31, "Ga", "gallium"},      {32, "Ge", "germanium"},
        {33, "As", "arsenic"},      {34, "Se", "selenium"},
        {35, "Br", "bromine"},      {36, "Kr", "krypton"},
        {37, "Rb", "rubidium"},     {38, "Sr", "strontium"},
        {39, "Y", "yttrium"},       {40, "Zr", "zirconium"},
        {41, "Nb", "niobium"},      {42, "Mo", "molybdenum"},
        {43, "Tc", "technetium"},   {44, "Ru", "ruthenium"},
        {45, "Rh", "rhodium"},      {46, "Pd", "palladium"},
        {47, "Ag", "silver"},       {48, "Cd", "cadmium"},
        {49, "In", "indium"},       {50, "Sn", "tin"},
        {51, "Sb", "antimony"},     {52, "Te", "tellurium"},
        {53, "I", "iodine"},        {54, "Xe", "xenon"},
        {55, "Cs", "caesium"},      {56, "Ba", "barium"},
        {57, "La", "lanthanum"},    {58, "Ce", "cerium"},
        {59, "Pr", "praseodymium"}, {60, "Nd", "neodymium"},
        {61, "Pm", "promethium"},   {62, "Sm", "samarium"},
        {63, "Eu", "europium"},     {64, "Gd", "gadolinium"},
        {65, "Tb", "terbium"},      {66, "Dy", "dysprosium"},
        {67, "Ho", "holmium"},      {68, "Er", "erbium"},
        {69, "Tm", "thulium"},      {70, "Yb", "ytterbium"},
        {71, "Lu", "lutetium"},     {72, "Hf", "hafnium"},
        {73, "Ta", "tantalum"},     {74, "W", "tungsten"},
        {75, "Re", "rhenium"},      {76, "Os", "osmium"},
        {77, "Ir", "iridium"},      {78, "Pt", "platinum"},
        {79, "Au", "gold"},         {80, "Hg", "mercury"},
        {81, "Tl", "thallium"},     {82, "Pb", "lead"},
        {83, "Bi", "bismuth"},      {84, "Po", "polonium"},
        {85, "At", "astatine"},     {86, "Rn", "radon"},
        {87, "Fr", "francium"},     {88, "Ra", "radium"},
        {89, "Ac", "actinium"},     {90, "Th", "thorium"},
        {91, "Pa", "protactinium"}, {92, "U", "uranium"},
        {93, "Np", "neptunium"},    {94, "Pu", "plutonium"},
        {95, "Am", "americium"},    {96, "Cm", "curium"},
        {97, "Bk", "berkelium"},    {98, "Cf", "californium"},
        {99, "Es", "einsteinium"},  {100, "Fm", "fermium"},
        {101, "Md", "mendelevium"}, {102, "No", "nobelium"},
        {103, "Lr", "lawrencium"},  {104, "Rf", "rutherfordium"},
        {105, "Ha", "hahnium"},     {106, "Sg", "seaborgium"},
        {107, "Bh", "bohrium"},     {108, "Hs", "hassium"},
        {109, "Mt", "meitnerium"},  {110, "Ds", "darmstadtium"},
        {111, "Rg", "roentgenium"}, {112, "Cn", "copernicium"},
        {113, "Nh", "nihonium"},    {114, "Fl", "flerovium"},
        {115, "Mc", "moscovium"},   {116, "Lv", "livermorium"},
        {117, "Ts", "tennessine"},  {118, "Og", "oganesson"},
  }};

  return elements;
}

bool ignore_case_equal(const std::string& a, const std::string& b) {
  return a.size() == b.size() &&
         std::equal(std::begin(a), std::end(a), std::begin(b), [](char ca, char cb) {
           return std::tolower(ca) == std::tolower(cb);
         });
}

bool is_element_symbol(const std::string& symbol) {
  return std::any_of(elements().begin(), elements().end(), [&symbol](const Element& e) {
    return ignore_case_equal(symbol, e.symbol);
  });
}

const Element& Element::by_symbol(const std::string& symbol) {
  auto res = std::find_if(
        std::begin(elements()), std::end(elements()),
        [&symbol](const Element& e) { return ignore_case_equal(symbol, e.symbol); });

  assert_throw(res != std::end(elements()),
               ExcUnknownElement("Element symbol \"" + symbol + "\" not known."));
  return *res;
}

bool is_atomic_number(unsigned int atomic_number) {
  return 0 < atomic_number && atomic_number <= elements().size();
}

const Element& Element::by_atomic_number(unsigned int atomic_number) {
  assert_throw(is_atomic_number(atomic_number),
               ExcUnknownElement("Only know atomic numbers in range [1," +
                                 std::to_string(elements().size()) + "]."));
  const auto& e = elements()[atomic_number - 1];
  assert_dbg(e.atomic_number == atomic_number, krims::ExcInternalError());
  return e;
}

}  // namespace gint
