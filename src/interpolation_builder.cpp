#include "interpolation_builder.hpp"
#include "real_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct LinearInterpolationBuilder : public InterpolationBuilder
      {
        virtual ~LinearInterpolationBuilder( void ) = default;
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
        }
      };
      
      struct NaturalCubicSplineInterpolationBuilder : public InterpolationBuilder
      {
        virtual ~NaturalCubicSplineInterpolationBuilder( void ) = default;
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createNaturalCubicSplineInterpolatedFunction( xValues, yValues );
        }
      };
      
      struct PiecewiseConstantRightInterpolationBuilder : public InterpolationBuilder
      {
        virtual ~PiecewiseConstantRightInterpolationBuilder( void ) = default;
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction( xValues, yValues );
        }
      };
    }

    beagle::interp_builder_ptr_t
    InterpolationBuilder::linear( void )
    {
      return std::make_shared<impl::LinearInterpolationBuilder>();
    }

    beagle::interp_builder_ptr_t
    InterpolationBuilder::naturalCubicSpline( void )
    {
      return std::make_shared<impl::NaturalCubicSplineInterpolationBuilder>();
    }

    beagle::interp_builder_ptr_t
    InterpolationBuilder::piecewiseConstantRight( void )
    {
      return std::make_shared<impl::PiecewiseConstantRightInterpolationBuilder>();
    }
  }
}