#pragma once
#include "Integral.hh"
#include "config.hh"
#include <krims/ParameterMap.hh>

namespace gint {

typedef enum { REAL_ATOMIC, COMPLEX_ATOMIC, REAL_MOLECULAR, COMPLEX_MOLECULAR } OrbitalType;

template <OrbitalType otype>
class IntegralCollectionBase {
public:
  constexpr static OrbitalType orbital_type = otype;

  typedef typename std::conditional<otype == COMPLEX_MOLECULAR, complex_stored_mtx_type,
                                    real_stored_mtx_type>::type stored_mtx_type;
  typedef Integral<stored_mtx_type> integral_matrix_type;

  virtual integral_matrix_type lookup_integral(const std::string& integral_id) const = 0;
};

template <OrbitalType otype>
using create_collection_t =
      std::shared_ptr<IntegralCollectionBase<otype>>(const krims::ParameterMap& parameters);

}  // namespace gint
