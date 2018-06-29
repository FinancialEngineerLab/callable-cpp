#include "interpolation_builder.hpp"
#include "real_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct LinearWithFlatExtrapolationInterpolationBuilder : public InterpolationBuilder
      {
        LinearWithFlatExtrapolationInterpolationBuilder( void )
        { }
        virtual ~LinearWithFlatExtrapolationInterpolationBuilder( void )
        { }
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
        }
      };
      
      struct NaturalCubicSplineWithFlatExtrapolationInterpolationBuilder : public InterpolationBuilder
      {
        NaturalCubicSplineWithFlatExtrapolationInterpolationBuilder( void )
        { }
        virtual ~NaturalCubicSplineWithFlatExtrapolationInterpolationBuilder( void )
        { }
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createNaturalCubicSplineWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
        }
      };
    }

    InterpolationBuilder::InterpolationBuilder( void )
    { }

    InterpolationBuilder::~InterpolationBuilder( void )
    { }

    beagle::interp_builder_ptr_t
    InterpolationBuilder::linearWithFlatExtrapolation( void )
    {
      return std::make_shared<impl::LinearWithFlatExtrapolationInterpolationBuilder>();
    }

    beagle::interp_builder_ptr_t
    InterpolationBuilder::naturalCubicSplineWithFlatExtrapolation( void )
    {
      return std::make_shared<impl::NaturalCubicSplineWithFlatExtrapolationInterpolationBuilder>();
    }
  }
}