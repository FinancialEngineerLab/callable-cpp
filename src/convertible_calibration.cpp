#include "convertible_calibration.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"

namespace beagle
{
  namespace calibration
  {
    namespace util
    {
      namespace
      {
        void formModelParameters(double sigma,
                                 double intensity,
                                 double spot,
                                 double exponent,
                                 beagle::real_2d_function_ptr_t& drift,
                                 beagle::real_2d_function_ptr_t& volatility,
                                 beagle::real_2d_function_ptr_t& rate)
        {
          drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                          [=](double time, double price){ return intensity * std::pow(price / spot, -exponent); } );
          volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);
          rate = drift;
        }
      }

      beagle::andersen_buffum_param_t
      createCalibratedAndersenBuffumParameters(const beagle::real_function_ptr_t& forward,
                                               const beagle::real_function_ptr_t& discounting,
                                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                               double exponent,
                                               const beagle::dbl_vec_t& expiries,
                                               const beagle::andersen_buffum_param_t& quotes)
      {
        if (expiries.size() != quotes.size())
          throw(std::string("The number of expiries is not the same as the number of calibrations for Andersen-Buffum!"));

        // Set up the spatial finite difference grid and the initial condition
        beagle::real_2d_function_ptr_t drift;
        beagle::real_2d_function_ptr_t volatility;
        beagle::real_2d_function_ptr_t rate;
        formModelParameters(quotes.back().first,
                            quotes.back().second,
                            forward->value(0.0),
                            exponent,
                            drift,
                            volatility,
                            rate);

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

        // Now perform calibration by bootstrapping
        beagle::andersen_buffum_param_t params;
        for (double expiry : expiries)
        {
          
        }


        return beagle::andersen_buffum_param_t{};
      }
    }





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

      return beagle::real_function_ptr_coll_t{beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction(expiries, a),
                                              beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction(expiries, b)};
    }
  }
}