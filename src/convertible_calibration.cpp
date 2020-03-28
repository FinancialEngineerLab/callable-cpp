#include "real_function.hpp"

namespace beagle
{
  namespace math
  {
    beagle::real_function_ptr_coll_t
    RealFunction::createCalibratedAndersenBuffumParameters(const beagle::real_function_ptr_t& forward,
                                                           const beagle::real_function_ptr_t& discounting,
                                                           const beagle::real_2d_function_ptr_t& drift,
                                                           const beagle::real_2d_function_ptr_t& volatility,
                                                           const beagle::real_2d_function_ptr_t& rate,
                                                           const beagle::real_2d_function_ptr_t& recovery,
                                                           const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                                           const beagle::dbl_vec_t& expiries,
                                                           const beagle::two_dbl_t& quotes )
    {
      return beagle::real_function_ptr_coll_t{};
    }
  }
}