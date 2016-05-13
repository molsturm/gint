# Installs the cmake apckage information this project provides
#
# Requires the variable PackageModuleLocation to be set.

# Write a basic version file for gint
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${gint_BINARY_DIR}/gintConfigVersion.cmake"
	COMPATIBILITY AnyNewerVersion
)

# Adjust a configure file
configure_file(cmake/gintConfig.cmake.in
	"${gint_BINARY_DIR}/gintConfig.cmake"
	COPYONLY
)

# Set an export location:
install(FILES
	"${gint_BINARY_DIR}/gintConfig.cmake"
	"${gint_BINARY_DIR}/gintConfigVersion.cmake"
	DESTINATION "${PackageModuleLocation}/gint"
	COMPONENT devel
)

