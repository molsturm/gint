#
# The sources for the test executable
#
set(GINT_TESTS_SOURCES
  QuickAtomicCsDummyTest.cc
  QuickAtomicStatic14Test.cc
  QuickAtomicCsNaiveTest.cc
  QuickLibintTests.cc
  main.cc
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

# Add the test executables
FOREACH(target ${gint_tests_TARGETS})
	add_test(NAME ${target} COMMAND ${target})
ENDFOREACH()

