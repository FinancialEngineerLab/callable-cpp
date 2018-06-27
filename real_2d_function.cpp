#include "real_2d_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct TwoDimConstantFunction : public RealTwoDimFunction
      {
        TwoDimConstantFunction( double constant ) :
          m_Const( constant )
        { }
        virtual ~TwoDimConstantFunction( void )
        { }
      public:
        virtual double value( double argX,
                              double argY ) const override
        {
          return m_Const;
        }
      private:
        double m_Const;
      };
    }


    RealTwoDimFunction::RealTwoDimFunction( void )
    { }

    RealTwoDimFunction::~RealTwoDimFunction( void )
    { }

    beagle::real_2d_function_ptr_t
    RealTwoDimFunction::createTwoDimConstantFunction( double constant )
    {

      return std::make_shared<impl::TwoDimConstantFunction>( constant );
    }
  }
}