#ifndef REAL_2D_FUNCTION_HPP
#define REAL_2D_FUNCTION_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct RealTwoDimFunction
    {
      RealTwoDimFunction( void );
      virtual ~RealTwoDimFunction( void );
    public:
      virtual double value( double argX,
                            double argY ) const = 0;
    public:
      static beagle::real_2d_function_ptr_t createTwoDimConstantFunction( double constant );
      static beagle::real_2d_function_ptr_t createPiecewiseConstantRightFunction(
                                  const beagle::dbl_vec_t& params,
                                  const beagle::real_function_ptr_coll_t& funcs );
      static beagle::real_2d_function_ptr_t createBootstrappedLocalVolatilityFunction(
                                  const beagle::dbl_vec_t& expiries,
                                  const beagle::dbl_vec_t& initialGuesses,
                                  const beagle::dbl_vec_vec_t& strikesColl,
                                  const beagle::dbl_vec_vec_t& pricesColl,
                                  const beagle::pricer_ptr_t& forwardPricer,
                                  const beagle::payoff_ptr_t& payoff,
                                  const beagle::interp_builder_ptr_t& interp);
    };
  }
}


#endif