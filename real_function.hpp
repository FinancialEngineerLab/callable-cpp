#ifndef REAL_FUNCTION_HPP
#define REAL_FUNCTION_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct RealFunction
    {
      RealFunction( void );
      virtual ~RealFunction( void );
    public:
      virtual double value( double arg ) const = 0;
    public:
      static beagle::real_function_ptr_t createConstantFunction( double constant );
      static beagle::real_function_ptr_t createLinearWithFlatExtrapolationInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues );
      static beagle::real_function_ptr_t createNaturalCubicSplineWithFlatExtrapolationInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues );
    };
  }
}


#endif