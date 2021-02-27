#pragma once

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct InterpolationBuilder
    {
      virtual ~InterpolationBuilder( void ) = default;
    public:
      virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                        const dbl_vec_t& yValues ) const = 0;
    public:
      static beagle::interp_builder_ptr_t linear( void );
      static beagle::interp_builder_ptr_t naturalCubicSpline( void );
      static beagle::interp_builder_ptr_t piecewiseConstantRight( void );
    };
  }
}
