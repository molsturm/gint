#pragma once
#include "gint/config.hh"
#include <string>

namespace gint {

struct version {
  static int constexpr major{detail::version_major};
  static int constexpr minor{detail::version_minor};
  static int constexpr patch{detail::version_patch};

  // Return the version as a string
  static std::string version_string();
};

}  // namespace gint
