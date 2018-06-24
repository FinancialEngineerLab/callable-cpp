#ifndef UTIL_HPP
#define UTIL_HPP

#include <cmath>

namespace beagle
{
  namespace util
  {
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
  }
}

#endif