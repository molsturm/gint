## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the gint authors
##
## This file is part of gint.
##
## gint is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gint is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with gint. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

set(GINT_TESTS_SOURCES
  ReadBasisSetTest.cc
  BasisSetTests.cc
  BasisTests.cc
  IntegralLookupTests.cc

  SturmianTestData.cc
  CoulombSturmianSmallTest.cc

  GaussianTestData.cc
  GaussianSmallTest.cc

  main.cc

  # TODO Make a proper test out of this
  CoulombSturmianCoefficientTests.cc
)

drb_add_executable(gint_tests
	DBGSUFFIX  "_DEBUG"
	RELSUFFIX  "_RELEASE" 
	FILES ${GINT_TESTS_SOURCES}
)

drb_target_compile_definitions(ALL     gint_tests
	PRIVATE ${GINT_DEFINITIONS})
drb_target_compile_definitions(DEBUG   gint_tests
	PRIVATE ${GINT_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE gint_tests
	PRIVATE ${GINT_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     gint_tests
	${GINT_DEPENDENCIES} ${GINT_DEPENDENCIES_TEST})
drb_target_link_libraries(DEBUG   gint_tests
	${gint_DEBUG_TARGET} ${GINT_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE gint_tests
	${gint_RELEASE_TARGET} ${GINT_DEPENDENCIES_RELEASE})

drb_add_dependencies(ALL  gint_tests  ${DATA_DOWN_TARGET})

# Add the test executables
FOREACH(target ${gint_tests_TARGETS})
	add_test(NAME ${target} COMMAND ${target})
ENDFOREACH()

