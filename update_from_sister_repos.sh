#!/bin/bash

. update_from_sisters.lib.sh || exit 1

update_file "krims" ".clang-format"

update_file "krims" "external/get.lib.sh" || exit 1
update_file "krims" "external/get_rapidcheck.sh" || exit 1
update_file "linalgwrap" "external/get_krims.sh" || exit 1
update_file "gscf" "external/get_linalgwrap.sh" || exit 1

update_file "krims" "cmake/findRapidcheck.cmake" || exit 1
update_file "krims" "cmake/findCatch.cmake" || exit 1
update_file "linalgwrap" "cmake/findKrims.cmake" || exit 1
update_file "gscf" "cmake/findLinalgwrap.cmake" || exit 1
update_file "linalgwrap" "cmake/IncludeKrimsCmakeModule.cmake" || exit 1

update_file "krims" "templates/cc.template" "keep_header" || exit 1
update_file "krims" "templates/cmake.template" "keep_header" || exit 1
update_file "krims" "templates/hh.template" "keep_header" || exit 1
update_file "krims" "templates/README.md" || exit 1

update_file "linalgwrap" "update_from_sisters.lib.sh" || exit 1
