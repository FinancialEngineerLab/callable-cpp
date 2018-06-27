#include "interpolation_builder.hpp"

namespace beagle
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
        return beagle::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction( xValues, yValues );
      }
    };
  }

  InterpolationBuilder::InterpolationBuilder( void )
  { }

  InterpolationBuilder::~InterpolationBuilder( void )
  { }

  beagle::real_function_ptr_t
  InterpolationBuilder::linearWithFlatExtrapolation( void )
  {
    return std::make_shared<impl::LinearWithFlatExtrapolationInterpolationBuilder>();
  }
}