#pragma once

#include "fwd_decl.hpp"

namespace beagle
{
  namespace util
  {
    const double& pi( void );

    const double& rootTwo( void );

    const double& epsilon( void );

    double standardNormal( double arg );

    double cumulativeStandardNormal( double arg );

    double bsCall( double strike,
                   double forward,
                   double expiry,
                   double vol );

    double bsVega( double strike,
                   double forward,
                   double expiry,
                   double vol );

    double impliedBlackVolatility( double price,
                                   double strike,
                                   double forward,
                                   double expiry,
                                   const beagle::real_function_ptr_t& discounting );

    void tridiagonalSolve( beagle::dbl_vec_t& rhs,
                           beagle::dbl_vec_t& diag,
                           beagle::dbl_vec_t& upper,
                           beagle::dbl_vec_t& lower );

    void inverseMatrix( beagle::dbl_vec_vec_t& inputMat );
  }
}
