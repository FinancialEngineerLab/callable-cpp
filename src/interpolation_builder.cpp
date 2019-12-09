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
        LinearInterpolationBuilder( void )
        { }
        virtual ~LinearInterpolationBuilder( void )
        { }
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
        }
      };
      
      struct NaturalCubicSplineInterpolationBuilder : public InterpolationBuilder
      {
        NaturalCubicSplineInterpolationBuilder( void )
        { }
        virtual ~NaturalCubicSplineInterpolationBuilder( void )
        { }
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createNaturalCubicSplineWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
        }
      };
      
      struct PiecewiseConstantRightInterpolationBuilder : public InterpolationBuilder
      {
        PiecewiseConstantRightInterpolationBuilder( void )
        { }
        virtual ~PiecewiseConstantRightInterpolationBuilder( void )
        { }
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction( xValues, yValues );
        }
      };
    }

    InterpolationBuilder::InterpolationBuilder( void )
    { }

    InterpolationBuilder::~InterpolationBuilder( void )
    { }

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