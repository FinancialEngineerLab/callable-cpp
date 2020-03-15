#ifndef UTIL_HPP
#define UTIL_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace util
  {
    beagle::dbl_t pi( void );
    
    beagle::dbl_t rootTwo( void );

    beagle::dbl_t standardNormal( beagle::dbl_t arg );

    beagle::dbl_t cumulativeStandardNormal( beagle::dbl_t arg );

    beagle::dbl_t bsCall( beagle::dbl_t strike,
                   beagle::dbl_t forward,
                   beagle::dbl_t expiry,
                   beagle::dbl_t vol );

    beagle::dbl_t bsVega( beagle::dbl_t strike,
                   beagle::dbl_t forward,
                   beagle::dbl_t expiry,
                   beagle::dbl_t vol );

    // The indices for all the vectors run from 0 to N+1, with the zeroth and (N+1)-th 
    // elements used for boundary conditions, up to second order
    // d_0 * x_0         + u_0 * x_1     + l_0 * x_2         = y_0
    // u_{N+1} * x_{N-1} + l_{N+1} * x_N + d_{N+1} * x_{N+1} = y_{N+1}
    void tridiagonalSolve( beagle::dbl_vec_t& rhs,
                           beagle::dbl_vec_t& diag,
                           beagle::dbl_vec_t& upper,
                           beagle::dbl_vec_t& lower );

    void inverseMatrix( beagle::dbl_vec_vec_t& inputMat );
  }
}

#endif