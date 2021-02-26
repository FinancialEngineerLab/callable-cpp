#include "gtest/gtest.h"
#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "util.hpp"
#include "integration_method.hpp"

namespace beagle
{
  namespace test
  {
    TEST(test_exact_sabr, ClosedFormValuation)
    {
      double spot = 50;
      double r = .04;
      double q = .02;

      double a = .25;
      double b = .3;
      double c = -.3;
      double d = .45;

      double expiry = 1.;
      double strike = 50.;
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::put();
      beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption(expiry,
                                                                                               strike,
                                                                                               payoff);

      // put option valuation with discounting, funding, and volatility
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_function_ptr_t alpha = beagle::math::RealFunction::createConstantFunction(a);
      beagle::real_function_ptr_t beta = beagle::math::RealFunction::createConstantFunction(b);
      beagle::real_function_ptr_t rho = beagle::math::RealFunction::createConstantFunction(c);
      beagle::real_function_ptr_t nu = beagle::math::RealFunction::createConstantFunction(d);
      beagle::integration_method_ptr_t quad = beagle::math::IntegrationMethod::midPointIntegrationMethod(500);

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormExactSABREuropeanOptionPricer(forward,
                                                                                                           discounting,
                                                                                                           alpha,
                                                                                                           beta,
                                                                                                           rho,
                                                                                                           nu,
                                                                                                           true,
                                                                                                           quad);

      printf("\n************************************************************************************\n");
      printf("Valuation results for the exact SABR model:\n");
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
      printf("************************************************************************************\n\n");
    }
  }
}