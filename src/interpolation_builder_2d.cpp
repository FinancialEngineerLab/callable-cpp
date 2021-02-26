#include "interpolation_builder_2d.hpp"


namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct TwoOneDimInterpolationBuilder : public TwoDimInterpolationBuilder
      {
        TwoOneDimInterpolationBuilder( const beagle::interp_builder_ptr_t& interpX,
                                       const beagle::interp_builder_ptr_t& interpY ) :
          m_InterpX( interpX ),
          m_InterpY( interpY )
        { }
      public:
        virtual beagle::real_2d_function_ptr_t formTwoDimFunction() const override
        {
          throw(std::string("Cannot form a two dimensional function with two one dimensional interpolation builders!"));
          return real_2d_function_ptr_t();
        }
      private:
        beagle::interp_builder_ptr_t m_InterpX;
        beagle::interp_builder_ptr_t m_InterpY;
      };
    }

    beagle::interp_builder_2d_ptr_t
    TwoDimInterpolationBuilder::formTwoOneDimInterpolationBuilder(
                                                  const beagle::interp_builder_ptr_t& interpX,
                                                  const beagle::interp_builder_ptr_t& interpY )
    {
      return std::make_shared<impl::TwoOneDimInterpolationBuilder>( interpX, interpY );
    }
  }
}
