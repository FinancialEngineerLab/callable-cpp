#ifndef REAL_2D_FUNCTION_HPP
#define REAL_2D_FUNCTION_HPP

#include <memory>

namespace beagle
{

  struct RealTwoDimFunction;
  using real_2d_function_ptr_t = std::shared_ptr<RealTwoDimFunction>;

  struct RealTwoDimFunction
  {
    RealTwoDimFunction( void );
    virtual ~RealTwoDimFunction( void );
  public:
    virtual double value( double argX,
                          double argY ) const = 0;
  public:
    static real_2d_function_ptr_t createTwoDimConstantFunction( double constant );
  };
}


#endif