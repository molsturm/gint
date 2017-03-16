#pragma once
#include <complex>
#include <linalgwrap/Armadillo/ArmadilloMatrix.hh>
#include <linalgwrap/Base/Interfaces/MutableMemoryVector_i.hh>

namespace gint {

typedef double real_type;
typedef std::complex<real_type> complex_type;

typedef linalgwrap::ArmadilloMatrix<real_type> real_stored_mtx_type;
typedef linalgwrap::ArmadilloMatrix<complex_type> complex_stored_mtx_type;

// typedef linalgwrap::MultiVector<const linalgwrap::MutableMemoryVector_i<real_type> >
// const_real_multivector_type;
// typedef linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<real_type> >
// real_multivector_type;
typedef const real_stored_mtx_type const_real_multivector_type;
typedef real_stored_mtx_type real_multivector_type;

}  // namespace gint
