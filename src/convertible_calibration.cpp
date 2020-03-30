#include "real_function.hpp"

namespace beagle
{
  namespace math
  {





    // This function calibrates a model to a constant ATMF volatility and a constant risky spread.
    // See Andersen and Buffum, Calibration and Implementation of Convertible Bond Models, for details.
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
      beagle::dbl_vec_t a(expiries.size());
      beagle::dbl_vec_t b(expiries.size());

      // Set up the spatial finite difference grid and the initial condition
      beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                               forward,
                                                               discounting,
                                                               drift,
                                                               volatility,
                                                               rate,
                                                               settings);
      auto pODFP = dynamic_cast<beagle::valuation::mixins::OneDimFokkerPlanck*>(forwardPricer.get());

      beagle::dbl_vec_t stateVars;
      beagle::dbl_vec_t density;
      pODFP->formInitialCondition(expiries.back(), stateVars, density);

      return beagle::real_function_ptr_coll_t{beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction(expiries, a),
                                              beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction(expiries, b)};
    }
  }
}