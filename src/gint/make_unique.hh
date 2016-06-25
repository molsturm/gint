#pragma once
#include <memory>

namespace gint {
  
template <class T, class... Args>
std::unique_ptr<T> make_unique(Args&&... args){
  return std::unique_ptr<T>(new T(args...));
}

}
