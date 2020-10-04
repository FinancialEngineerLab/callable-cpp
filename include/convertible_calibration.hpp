#ifndef CONVERTIBLE_CALIBRATION_HPP
#define CONVERTIBLE_CALIBRATION_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace calibration
  {
    namespace util
    {
      // This function calibrates a model to a constant ATMF volatility and a constant risky spread.
      // See Andersen and Buffum, Calibration and Implementation of Convertible Bond Models, for details.
      beagle::andersen_buffum_curve_pair_t
      createCalibratedAndersenBuffumParameters(const beagle::real_function_ptr_t& forward,
                                               const beagle::real_function_ptr_t& discounting,
                                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                               double exponent,
                                               const beagle::dbl_vec_t& expiries,
                                               const beagle::andersen_buffum_param_t& quotes);

      // This function calibrates a model to an implied volatility surface together with credit spread
      beagle::andersen_buffum_curve_pair_t
      createCalibratedAndersenBuffumParameters(const beagle::real_function_ptr_t& forward,
                                               const beagle::real_function_ptr_t& discounting,
                                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                               double exponent,
                                               const beagle::volatility_smile_credit_spread_coll_t& quotes,
                                               const beagle::interp_builder_ptr_t& interp);
    }
  }
}

#endif