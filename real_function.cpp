#include "real_function.hpp"

namespace beagle
{

  namespace impl
  {
    struct ConstantFunction : public RealFunction
    {
      ConstantFunction( double constant ) :
        m_Const( constant )
      { }
      virtual ~ConstantFunction( void )
      { }
    public:
      virtual double value( double arg ) const override
      {
        return m_Const;
      }
    private:
      double m_Const;
    };
  }


  RealFunction::RealFunction( void )
  { }

  RealFunction::~RealFunction( void )
  { }

  real_function_ptr_t
  RealFunction::createConstantFunction( double constant )
  {

    return std::make_shared<impl::ConstantFunction>( constant );
  }
}