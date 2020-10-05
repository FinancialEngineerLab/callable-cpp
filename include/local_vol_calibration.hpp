#pragma once

#include "fwd_decl.hpp"

namespace beagle
{
  namespace calibration
  {
    namespace util
    {
      beagle::real_2d_function_ptr_t
      createCalibratedLocalVolatilitySurface(const beagle::real_function_ptr_t& forward,
                                             const beagle::real_function_ptr_t& discounting,
                                             const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                             const beagle::volatility_smile_coll_t& volSmiles);
    }
  }
}
