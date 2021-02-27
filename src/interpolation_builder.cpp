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
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
        }
      };

      struct NaturalCubicSplineInterpolationBuilder : public InterpolationBuilder
      {
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createNaturalCubicSplineWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
        }
      };

      struct PiecewiseConstantRightInterpolationBuilder : public InterpolationBuilder
      {
      public:
        virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                          const dbl_vec_t& yValues ) const override
        {
          return beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction( xValues, yValues );
        }
      };
    }

    const beagle::interp_builder_ptr_t&
    InterpolationBuilder::linear( void )
    {
      static beagle::interp_builder_ptr_t instance = std::make_shared<impl::LinearInterpolationBuilder>();
      return instance;
    }

    const beagle::interp_builder_ptr_t&
    InterpolationBuilder::naturalCubicSpline( void )
    {
      static beagle::interp_builder_ptr_t instance = std::make_shared<impl::NaturalCubicSplineInterpolationBuilder>();
      return instance;
    }

    const beagle::interp_builder_ptr_t&
    InterpolationBuilder::piecewiseConstantRight( void )
    {
      static beagle::interp_builder_ptr_t instance = std::make_shared<impl::PiecewiseConstantRightInterpolationBuilder>();
      return instance;
    }
  }
}