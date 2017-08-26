## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the gint authors
##
## This file is part of gint.
##
## gint is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gint is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with gint. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

set(GINT_SOURCES
  version.cc
  qnum.cc

  Atom.cc
  CoefficientContainer.cc
  Element.cc
  elements.cc
  ERITensor_i.cc
  Integral.cc
  IntegralCoreBase.cc
  IntegralLookup.cc
  IntegralLookupKeys.cc
  IntegralType.cc
  IntegralUpdateKeys.cc
  Structure.cc

  sturmian/atomic/nlm_order/ERICore.cc
  sturmian/atomic/nlm_order/ERICoreHighPrecision.cc
  sturmian/atomic/nlm_order/ERITensor.cc
  sturmian/atomic/nlm_order/IntegralCollectionBase.cc
  sturmian/atomic/nlm_order/SturmintSystem.cc
  sturmian/atomic/cs_dummy.cc
  sturmian/atomic/cs_naive.cc
  sturmian/atomic/cs_reference.cc
  sturmian/atomic/cs_reference_pc.cc
  sturmian/atomic/cs_static14.cc
  sturmian/atomic/NlmBasis.cc
  sturmian/atomic/Static14Data.cc

  gaussian/libcint.cc
  gaussian/libint.cc
  gaussian/Basis.cc
  gaussian/BasisSet.cc
  gaussian/read_basisset.cc
  gaussian/IntegralLookupKeys.cc
)

# Write the config file config.hh (into the binary dir)
configure_file("config.hh.in" "config.hh")

# Add the libraries for Debug and Release builds:
drb_add_library(gint
	DBGSUFFIX  ".g"
	FILES ${GINT_SOURCES}
)

drb_set_target_properties(
	ALL gint
	PROPERTIES
	VERSION "${PROJECT_VERSION}"
)

drb_target_link_libraries(ALL     gint ${GINT_DEPENDENCIES})
drb_target_link_libraries(DEBUG   gint ${GINT_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE gint ${GINT_DEPENDENCIES_RELEASE})

# set the installation property for these targets:
install(TARGETS ${gint_TARGETS}
	EXPORT GintTargets
	DESTINATION lib
	COMPONENT library
)

# install the includes:
install(DIRECTORY .
	DESTINATION "include/gint"
	COMPONENT devel
	FILES_MATCHING PATTERN "*.hh"
)
install(FILES "${gint_BINARY_DIR}/src/gint/config.hh"
	DESTINATION "include/gint"
	COMPONENT devel
)

# Export the target specifications for gint.
export(EXPORT GintTargets
	FILE "${gint_BINARY_DIR}/gintTargets.cmake"
	NAMESPACE Upstream::
)

# and install them appropriately:
install(EXPORT GintTargets
	FILE "gintTargets.cmake"
	NAMESPACE Upstream::
	DESTINATION "share/cmake/gint"
)

