#ifndef REAL_2D_FUNCTION_HPP
#define REAL_2D_FUNCTION_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct RealTwoDimFunction
    {
      virtual ~RealTwoDimFunction( void ) = default;
    public:
      virtual double value( double argX,
                            double argY ) const = 0;
    public:
      static beagle::real_2d_function_ptr_t createBinaryFunction( const beagle::real_2d_func_t& func );
      static beagle::real_2d_function_ptr_t createTwoDimConstantFunction( double constant );
      static beagle::real_2d_function_ptr_t createPiecewiseConstantRightFunction(
                                  const beagle::dbl_vec_t& params,
                                  const beagle::real_function_ptr_coll_t& funcs );
    };
  }
}


#endif