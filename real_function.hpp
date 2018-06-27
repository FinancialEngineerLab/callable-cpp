#ifndef REAL_FUNCTION_HPP
#define REAL_FUNCTION_HPP

#include <memory>
#include <vector>

namespace beagle
{

  struct RealFunction;
  using real_function_ptr_t = std::shared_ptr<RealFunction>;
  using dbl_vec_t = std::vector<double>;

  struct RealFunction
  {
    RealFunction( void );
    virtual ~RealFunction( void );
  public:
    virtual double value( double arg ) const = 0;
  public:
    static real_function_ptr_t createConstantFunction( double constant );
    static real_function_ptr_t createLinearWithFlatExtrapolationInterpolatedFunction( 
                                                       const dbl_vec_t& xValues,
                                                       const dbl_vec_t& yValues );
  };
}


#endif