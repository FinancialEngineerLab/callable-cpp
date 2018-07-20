#ifndef UTIL_HPP
#define UTIL_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace util
  {
    double pi( void );
    
    double rootTwo( void );

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

    void tridiagonalSolve( beagle::dbl_vec_t& rhs,
                           beagle::dbl_vec_t& diag,
                           beagle::dbl_vec_t& upper,
                           beagle::dbl_vec_t& lower );

    void inverseMatrix( beagle::dbl_vec_vec_t& inputMat );
  }
}

#endif