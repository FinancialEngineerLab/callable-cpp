#ifndef REAL_FUNCTION_HPP
#define REAL_FUNCTION_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct RealFunction
    {
      virtual ~RealFunction( void ) = default;
    public:
      virtual double value( double arg ) const = 0;
    public:
      static beagle::real_function_ptr_t createConstantFunction( double constant );
      static beagle::real_function_ptr_t createUnaryFunction( const beagle::real_func_t& func );
      static beagle::real_function_ptr_t createCompositeFunction( const beagle::real_function_ptr_t& f,
                                                                  const beagle::real_function_ptr_t& g );
      static beagle::real_function_ptr_t createLinearWithFlatExtrapolationInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues );
      static beagle::real_function_ptr_t createNaturalCubicSplineWithFlatExtrapolationInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues );
      static beagle::real_function_ptr_t createPiecewiseConstantRightInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues );
      static beagle::real_function_ptr_t createContinuousForwardAssetPriceFunction(
                                                         double spot,
                                                         const beagle::real_function_ptr_t& funding );
      static beagle::real_function_ptr_t createGeneralForwardAssetPriceFunction(
                                                         double spot,
                                                         const beagle::real_function_ptr_t& funding,
                                                         const beagle::dividend_schedule_t& dividends,
                                                         const beagle::dividend_policy_ptr_t& policy);
    };

    namespace mixins
    {
      struct InterpolationParameters
      {
        virtual ~InterpolationParameters( void ) = default;
      public:
        virtual const beagle::dbl_vec_t& xParameters( void ) const = 0;
        virtual const beagle::dbl_vec_t& yParameters( void ) const = 0;
      };

      struct ContainsDividends
      {
        virtual ~ContainsDividends( void ) = default;
      public:
        virtual const beagle::dividend_schedule_t& dividendSchedule( void ) const = 0;
        virtual const beagle::dividend_policy_ptr_t& dividendPolicy( void ) const = 0;
        virtual const beagle::cum_ex_dividend_prices_t& cumAndExDividendForwards( void ) const = 0;
      };
    }

    namespace impl
    {
      struct InterpolatedFunction : public RealFunction,
                                    public beagle::math::mixins::InterpolationParameters
      {
        InterpolatedFunction( const beagle::dbl_vec_t& xValues,
                              const beagle::dbl_vec_t& yValues );
        virtual ~InterpolatedFunction( void ) = default;
      public:
        virtual double value( double arg ) const = 0;
        virtual const beagle::dbl_vec_t& xParameters( void ) const override;
        virtual const beagle::dbl_vec_t& yParameters( void ) const override;
      private:
        beagle::dbl_vec_t m_XValues;
        beagle::dbl_vec_t m_YValues;
      };
    }
  }
}


#endif