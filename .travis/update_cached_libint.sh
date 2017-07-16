#!/bin/sh

# The directory where the libint cached files should live
CACHEDIR="$HOME/cache_libint"

# Place where the libint files are installed by ExternalProject
LIBINT_INSTALL_DIR="${TRAVIS_BUILD_DIR}/build/external/libint"

# -------------------------------------------------------------

# Replace cache by our build.
if [ -d "$LIBINT_INSTALL_DIR" ]; then
	echo "Updating libint cache at $CACHEDIR"

	rm -rf "$CACHEDIR"
	mkdir "$CACHEDIR"

	for subdir in include lib share; do
		cp -a $LIBINT_INSTALL_DIR/$subdir $CACHEDIR
	done
fi

