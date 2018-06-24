#ifndef REAL_FUNCTION_HPP
#define REAL_FUNCTION_HPP

#include <memory>

namespace beagle
{

  struct RealFunction;
  using real_function_ptr_t = std::shared_ptr<RealFunction>;

  struct RealFunction
  {
    RealFunction( void );
    virtual ~RealFunction( void );
  public:
    virtual double value( double arg ) const = 0;
  public:
    static real_function_ptr_t createConstantFunction( double constant );
  };
}


#endif