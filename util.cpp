#include "util.hpp"

namespace beagle
{
  namespace util
  {
    const double PI = 3.14159265358979323846;
    const double rootTwo = 1.41421356237;

    double standardNormal( double arg )
    {
      return std::exp( -.5*arg*arg ) / std::sqrt( 2 * PI );
    }

    double cumulativeStandardNormal( double arg )
    {
      return .5 * ( 1. + std::erf( arg / rootTwo ) );
    }

    double bsCall( double strike,
                   double forward,
                   double expiry,
                   double vol )
    {
      double moneyness = std::log( forward / strike );
      double totalDev = vol * std::sqrt( expiry );
      double dOne = moneyness / totalDev + .5 * totalDev;
      double dTwo = dOne - totalDev;

      return forward * cumulativeStandardNormal( dOne ) - strike * cumulativeStandardNormal( dTwo );
    }

    double bsVega( double strike,
                   double forward,
                   double expiry,
                   double vol )
    {
      double moneyness = std::log( forward / strike );
      double totalDev = vol * std::sqrt( expiry );
      double dOne = moneyness / totalDev + .5 * totalDev;

      return forward * standardNormal( dOne ) * std::sqrt( expiry );
    }
  }
}