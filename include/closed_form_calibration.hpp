#pragma once

#include "fwd_decl.hpp"

namespace beagle
{
  namespace calibration
  {
    namespace util
    {
      beagle::real_function_ptr_coll_t
      createCalibratedClosedFormEuropeanOptionPricerParameters(const beagle::real_function_ptr_t& forward,
                                                               const beagle::real_function_ptr_t& discounting,
                                                               const beagle::volatility_smile_coll_t& volSmiles,
                                                               const beagle::pricer_ptr_t& pricer,
                                                               const beagle::dbl_vec_t& guesses,
                                                               const beagle::calibration_bound_constraint_coll_t& constraints,
                                                               const beagle::interp_builder_ptr_coll_t& interps,
                                                               double tolerance = 1e-12);
    }
  }
}
