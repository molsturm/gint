#pragma once
#include <memory>

namespace gint {
// TODO: This should live in krams

/** Function to mimic the c++14 function std::make_unique
 *
 * Constructs a std::unique_ptr
 **/
template <class T, class... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(args...));
}

// TODO default to std::make_unique if c++14 is avail
}
