#ifndef INTERPOLATION_BUILDER_2D_HPP
#define INTERPOLATION_BUILDER_2D_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct TwoDimInterpolationBuilder
    {
      virtual ~TwoDimInterpolationBuilder( void ) = default;
    public:
      virtual beagle::real_2d_function_ptr_t formTwoDimFunction() const = 0;
    public:
      static beagle::interp_builder_2d_ptr_t formTwoOneDimInterpolationBuilder(
                                                  const beagle::interp_builder_ptr_t& interpX,
                                                  const beagle::interp_builder_ptr_t& interpY );
    };
  }
}



#endif