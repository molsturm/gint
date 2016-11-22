#pragma once

#include <complex>
typedef double real_type;
typedef std::complex<real_type> complex_type;

#include <linalgwrap/Armadillo/ArmadilloMatrix.hh>
typedef linalgwrap::ArmadilloMatrix<real_type>    real_stored_mtx_type;
typedef linalgwrap::ArmadilloMatrix<complex_type> complex_stored_mtx_type;


