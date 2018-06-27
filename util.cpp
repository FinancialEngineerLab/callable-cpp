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

    void tridiagonalSolve( beagle::dbl_vec_t& rhs,
                           beagle::dbl_vec_t& diag,
                           beagle::dbl_vec_t& upper,
                           beagle::dbl_vec_t& lower )
    {
      int size = rhs.size();
      for (int i = 0; i < size - 1; ++i)
      {
        diag[i+1] -= lower[i+1] / diag[i] * upper[i];
        rhs[i+1] -= lower[i+1] / diag[i] * rhs[i];
      }

      rhs[size-1] /= diag[size-1];

      for (int i = 0; i < size - 1; ++i)
      {
        int j = size - 2 - i;
        rhs[j] = (rhs[j] - upper[j] * rhs[j+1]) / diag[j];
      }
    }
  }
}