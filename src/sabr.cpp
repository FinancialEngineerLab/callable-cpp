#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "util.hpp"
#include "sabr.hpp"

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
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
      beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption(expiry,
                                                                                               strike,
                                                                                               payoff);

      // Call option valuation with discounting, funding, and volatility
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

      beagle::pricer_ptr_t sabr = beagle::valuation::Pricer::formClosedFormSABREuropeanOptionPricer(forward,
                                                                                                    discounting,
                                                                                                    alpha,
                                                                                                    beta,
                                                                                                    rho,
                                                                                                    nu);
      printf("%.9f\n", sabr->value(euroOption));
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
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
      beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption(expiry,
                                                                                               strike,
                                                                                               payoff);

      // Call option valuation with discounting, funding, and volatility
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

      beagle::pricer_ptr_t sabr = beagle::valuation::Pricer::formClosedFormShiftedSABREuropeanOptionPricer(forward,
                                                                                                           discounting,
                                                                                                           alpha,
                                                                                                           beta,
                                                                                                           rho,
                                                                                                           nu,
                                                                                                           shift);
      printf("%.9f\n", sabr->value(euroOption));
    }
  }
}