#include "gtest/gtest.h"
#include "pricer.hpp"
#include "option.hpp"
#include "bond.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"

namespace beagle
{
  namespace test
  {
    TEST(test_cev, ClosedFormValuation)
    {
      beagle::discrete_dividend_schedule_t dividends;
      // dividends.emplace_back( 0.5, 6.0 );
      // dividends.emplace_back( 1.5, 6.5 );
      // dividends.emplace_back( 2.5, 7.0 );
      // dividends.emplace_back( 3.5, 7.5 );
      // dividends.emplace_back( 4.5, 8.0 );
      // dividends.emplace_back( 5.5, 8.0 );
      // dividends.emplace_back( 6.5, 8.0 );

      double expiry = 1.;
      double strike = 100.;
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::put();
      beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                strike,
                                                                                                payoff );
      beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
                                                                                                strike,
                                                                                                payoff );

      double spot = 100.;
      double rate = .01;

      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return 1;});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot, discounting);

      // CEV modeling
      double alpha = .3;
      double beta = .5;
      beagle::real_function_ptr_t localVolFunction = beagle::math::RealFunction::createUnaryFunction(
                                                                  [alpha, beta](double x) { return alpha * std::pow(x, beta - 1.); } );
      beagle::real_2d_function_ptr_t localVolSurface =
                beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( beagle::dbl_vec_t(1U, 1.),
                                                                                        beagle::real_function_ptr_coll_t(1U, localVolFunction) );

      try
      {
        beagle::real_2d_function_ptr_t cev = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return alpha * std::pow(price, beta - 1.); } );
        beagle::pricer_ptr_t odbpop2  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                  forward,
                                                                  discounting,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  cev,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(500, 1001, 7.5) );
        beagle::pricer_ptr_t odfpeop2 = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
                                                                  forward,
                                                                  discounting,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  cev,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(500, 1001, 7.5) );

        printf("\n************************************************************************************\n");
        std::cout << "European option price (FD-B) is: " << odbpop2->value( euroOption ) << std::endl;
        std::cout << "American option price (FD-B) is: " << odbpop2->value( amerOption ) << std::endl;
        std::cout << "European option price (FD-F) is: " << odfpeop2->value( euroOption ) << std::endl;
        printf("************************************************************************************\n\n");
      }
      catch (const std::string& what)
      {
        std::cout << "A valuation error has occurred -- " << what << std::endl;
      }
    }
  }
}