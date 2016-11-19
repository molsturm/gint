#include "gint/version.hh"
#include <sstream>

namespace gint {

std::string version::version_string() {
  std::stringstream ss;
  ss << major << "." << minor << "." << patch;
  return ss.str();
}

}  // namespace gint
