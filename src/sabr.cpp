#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "util.hpp"
#include "sabr.hpp"
#include "integration_method.hpp"

#include <iostream>

namespace beagle
{
  namespace test
  {
    void test_sabr( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;

      double a = .25;
      double b = .5;
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

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormSABREuropeanOptionPricer(forward,
                                                                                                      discounting,
                                                                                                      alpha,
                                                                                                      beta,
                                                                                                      rho,
                                                                                                      nu);
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }

    void test_shifted_sabr( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;

      double a = .25;
      double b = .5;
      double c = -.3;
      double d = .45;
      double shift = 1.;

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

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormShiftedSABREuropeanOptionPricer(forward,
                                                                                                           discounting,
                                                                                                           alpha,
                                                                                                           beta,
                                                                                                           rho,
                                                                                                           nu,
                                                                                                           shift);
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }

    void test_cev(void)
    {
      double spot = 50.;
      double r = .04;
      double q = .02;

      double b = .25;
      double d = .3;

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
      beagle::real_function_ptr_t beta = beagle::math::RealFunction::createConstantFunction(b);
      beagle::real_function_ptr_t nu = beagle::math::RealFunction::createConstantFunction(d);
      beagle::integration_method_ptr_t quad = beagle::math::IntegrationMethod::midPointIntegrationMethod(500);

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormExactCEVEuropeanOptionPricer(forward,
                                                                                                       discounting,
                                                                                                       beta,
                                                                                                       nu,
                                                                                                       quad);
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }

    void test_free_boundary_cev( void )
    {
      double spot = 50.;
      double r = .04;
      double q = .02;

      double b = .25;
      double d = .3;

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
      beagle::real_function_ptr_t beta = beagle::math::RealFunction::createConstantFunction(b);
      beagle::real_function_ptr_t nu = beagle::math::RealFunction::createConstantFunction(d);
      beagle::integration_method_ptr_t quad = beagle::math::IntegrationMethod::midPointIntegrationMethod(500);

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormFreeBoundaryCEVEuropeanOptionPricer(forward,
                                                                                                              discounting,
                                                                                                              beta,
                                                                                                              nu,
                                                                                                              quad);
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }

    void test_exact_sabr(void)
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
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }

    void test_free_boundary_sabr(void)
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

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormFreeBoundarySABREuropeanOptionPricer(forward,
                                                                                                                discounting,
                                                                                                                alpha,
                                                                                                                beta,
                                                                                                                rho,
                                                                                                                nu,
                                                                                                                true,
                                                                                                                quad);
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }

    void test_normal_free_boundary_sabr(void)
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

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormNormalFreeBoundarySABREuropeanOptionPricer(forward,
                                                                                                                      discounting,
                                                                                                                      alpha,
                                                                                                                      rho,
                                                                                                                      nu,
                                                                                                                      true,
                                                                                                                      quad);
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }

    void test_normal_enhanced_free_boundary_sabr(void)
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

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer(forward,
                                                                                                                              discounting,
                                                                                                                              alpha,
                                                                                                                              beta,
                                                                                                                              rho,
                                                                                                                              nu,
                                                                                                                              true,
                                                                                                                              quad);
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
    }
  }
}