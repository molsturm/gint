#
# set the source files.
#
set(GINT_SOURCES
  version.cc

  IntegralUpdateKeys.cc
  Integral.cc
  IntegralCoreBase.cc
  IntegralLookup.cc
  IntegralType.cc
  CoefficientContainer.cc

  Atom.cc
  Element.cc
  Structure.cc

  sturmian/atomic/nlm_order/NlmBasis.cc
  sturmian/atomic/nlm_order/ERICore.cc
  sturmian/atomic/cs_dummy.cc
  sturmian/atomic/cs_naive.cc
  sturmian/atomic/cs_reference.cc
  sturmian/atomic/cs_reference_pc.cc
  sturmian/atomic/cs_static14.cc
  sturmian/atomic/Static14Data.cc

  gaussian/libint.cc
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

drb_target_compile_definitions(ALL     gint PUBLIC ${GINT_DEFINITIONS})
drb_target_compile_definitions(DEBUG   gint PUBLIC ${GINT_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE gint PUBLIC ${GINT_DEFINITIONS_RELEASE})

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
	DESTINATION "${PackageModuleLocation}/gint"
)

