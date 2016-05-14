#pragma once
#include <linalgwrap/Subscribable.hh>

namespace gint {

template <typename StoredMatrix>
class IntegralCoreBase : public linalgwrap::Subscribable {
public:
  typedef StoredMatrix stored_matrix_type;

  /** Friendly name of the integral */
  std::string name() const { return "make me pure"; }

  /** Identifier of the integral */
  std::string id() const { return "make me pure"; }

  // TODO specify interface and implement
};

}  // namespace gint
