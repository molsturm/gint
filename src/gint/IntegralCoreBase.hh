#pragma once
#include <linalgwrap/Subscribable.hh>

namespace gint {

template <typename StoredMatrix>
class IntegralCoreBase : public linalgwrap::Subscribable {
public:
  typedef StoredMatrix stored_matrix_type;

  // TODO specify interface and implement
};

}  // namespace gint
