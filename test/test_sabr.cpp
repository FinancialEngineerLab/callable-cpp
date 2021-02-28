#include "gtest/gtest.h"
#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "util.hpp"
#include "integration_method.hpp"
#include "closed_form_calibration.hpp"
#include "calibration_constraint.hpp"
#include "interpolation_builder.hpp"

namespace beagle
{
  namespace test
  {
    TEST(test_sabr, ClosedFormValuation)
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

      printf("\n************************************************************************************\n");
      printf("Valuation results for the SABR model:\n");
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
      printf("************************************************************************************\n\n");
    }

    TEST(test_sabr, Calibration)
    {
      double spot = 50;
      double r = .04;
      double q = .02;

      double a = .25;
      double b = .5;
      double c = -.3;
      double d = .45;

      // Implied volatilities are generated with:
      // alpha =  .25
      // beta  =  .3
      // rho   = -.3
      // nu    =  .45
      beagle::volatility_smile_coll_t smiles;
      double expiry = 1.;
      beagle::dbl_vec_t strikes = {48., 49., 50., 51., 52.};
      beagle::dbl_vec_t vols = {0.023776422, 0.020896676, 0.01822638, 0.016188335, 0.01553533};
      smiles.emplace_back(expiry, std::make_pair(strikes, vols));

      expiry = 2.;
      strikes = {47., 48., 49., 50., 51., 52., 53.};
      vols = {0.026720955, 0.023776422, 0.020896676, 0.01822638, 0.016188335, 0.01553533, 0.016359933};
      smiles.emplace_back(expiry, std::make_pair(strikes, vols));

      // put option valuation with discounting, funding, and volatility
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));

      beagle::dbl_vec_t guesses = {.2, .4, -.2, .4};
      beagle::real_function_ptr_t alpha = beagle::math::RealFunction::createConstantFunction(guesses[0]);
      beagle::real_function_ptr_t beta = beagle::math::RealFunction::createConstantFunction(guesses[1]);
      beagle::real_function_ptr_t rho = beagle::math::RealFunction::createConstantFunction(guesses[2]);
      beagle::real_function_ptr_t nu = beagle::math::RealFunction::createConstantFunction(guesses[3]);

      beagle::pricer_ptr_t pricer = beagle::valuation::Pricer::formClosedFormSABREuropeanOptionPricer(forward,
                                                                                                      discounting,
                                                                                                      alpha,
                                                                                                      beta,
                                                                                                      rho,
                                                                                                      nu);

      beagle::calibration_bound_constraint_coll_t constraints;
      constraints.push_back( beagle::calibration::CalibrationBoundConstraint::lowerBoundCalibrationConstraint(.0001) );
      constraints.push_back( beagle::calibration::CalibrationBoundConstraint::twoSidedBoundCalibrationConstraint(.0001, .9999) );
      constraints.push_back( beagle::calibration::CalibrationBoundConstraint::twoSidedBoundCalibrationConstraint(-.9999, .9999) );
      constraints.push_back( beagle::calibration::CalibrationBoundConstraint::lowerBoundCalibrationConstraint(.0001) );

      beagle::interp_builder_ptr_coll_t interps(4, beagle::math::InterpolationBuilder::linear());

      beagle::real_function_ptr_coll_t results = beagle::calibration::util::createCalibratedClosedFormEuropeanOptionPricerParameters(forward,
                                                                                                                                     discounting,
                                                                                                                                     smiles,
                                                                                                                                     pricer,
                                                                                                                                     guesses,
                                                                                                                                     constraints,
                                                                                                                                     interps);

      printf("\n************************************************************************************\n");
      auto pCF = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      beagle::pricer_ptr_t calibPricer = pCF->updateModelParameters(results);
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::put();
      printf("\n\nThe calibration result is:\n");
      for (int i = 0; i < smiles.size(); ++i)
      {
        double expiry = smiles[i].first;
        const beagle::dbl_vec_t& strikes = smiles[i].second.first;
        const beagle::dbl_vec_t& vols = smiles[i].second.second;
        for (int j = 0; j < strikes.size(); ++j)
        {
          beagle::product_ptr_t option = beagle::product::option::Option::createEuropeanOption(expiry, strikes[j], payoff);
          double target = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(forward, discounting, vols[j])->value(option);
          double source = calibPricer->value(option);
          printf("%.9f  %.9f\n", target, source);
        }
        printf("\n");
      }
      printf("End of showing calibration results.\n\n\n");
      printf("************************************************************************************\n\n");
    }

    TEST(test_shifted_sabr, ClosedFormValuation)
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

      printf("\n************************************************************************************\n");
      printf("Valuation results for the shifted SABR model:\n");
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
      printf("************************************************************************************\n\n");
    }

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

    TEST(test_free_boundary_sabr, ClosedFormValuation)
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

      printf("\n************************************************************************************\n");
      printf("Valuation results for the free boundary SABR model:\n");
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
      printf("************************************************************************************\n\n");
    }

    TEST(test_normal_enhanced_free_boundary_sabr, ClosedFormValuation)
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

      printf("\n************************************************************************************\n");
      printf("Valuation results for the normal enhanced free boundary SABR model:\n");
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
      printf("************************************************************************************\n\n");
    }

    TEST(test_normal_free_boundary_sabr, ClosedFormValuation)
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

      printf("\n************************************************************************************\n");
      printf("Valuation results for the normal free boundary SABR model:\n");
      printf("The option price is:          %.9f\n", pricer->value(euroOption));
      auto pIV = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
      printf("The implied BS volatility is: %.9f\n", pIV->impliedBlackScholesVolatility(euroOption));
      printf("************************************************************************************\n\n");
    }
  }
}