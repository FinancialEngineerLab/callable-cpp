#ifndef INTERPOLATION_BUILDER_HPP
#define INTERPOLATION_BUILDER_HPP

#include <memory>
#include <vector>
#include "real_function.hpp"

namespace beagle
{
  struct InterpolationBuilder;
  using interp_builder_ptr_t = std::shared_ptr<InterpolationBuilder>;
  using dbl_vec_t = std::vector<double>;

  struct InterpolationBuilder
  {
    InterpolationBuilder( void );
    virtual ~InterpolationBuilder( void );
  public:
    virtual beagle::real_function_ptr_t formFunction( const dbl_vec_t& xValues,
                                                      const dbl_vec_t& yValues ) const = 0;
  public:
    static interp_builder_ptr_t linearWithFlatExtrapolation( void );
  };
}



#endif