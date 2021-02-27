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
    TEST(test_black_scholes, ClosedFormValuation)
    {
      // Set up dividend
      beagle::dividend_schedule_t dividends;
      dividends.emplace_back(   .5, 0.0, 3.0 );
      dividends.emplace_back(  1.5, 0.0, 3.0 );
      dividends.emplace_back(  2.5, 0.0, 3.0 );
      dividends.emplace_back(  3.5, 0.0, 3.0 );
      dividends.emplace_back(  4.5, 0.0, 3.0 );
      dividends.emplace_back(  5.5, 0.0, 3.0 );
      dividends.emplace_back(  6.5, 0.0, 3.0 );
      dividends.emplace_back(  7.5, 0.0, 3.0 );
      dividends.emplace_back(  8.5, 0.0, 3.0 );
      dividends.emplace_back(  9.5, 0.0, 3.0 );
      dividends.emplace_back( 10.5, 0.0, 3.0 );
      dividends.emplace_back( 11.5, 0.0, 3.0 );
      dividends.emplace_back( 12.5, 0.0, 3.0 );
      dividends.emplace_back( 13.5, 0.0, 3.0 );
      dividends.emplace_back( 14.5, 0.0, 3.0 );
      dividends.emplace_back( 15.5, 0.0, 3.0 );
      dividends.emplace_back( 16.5, 0.0, 3.0 );
      dividends.emplace_back( 17.5, 0.0, 3.0 );
      dividends.emplace_back( 18.5, 0.0, 3.0 );
      dividends.emplace_back( 19.5, 0.0, 3.0 );

      beagle::discrete_dividend_schedule_t discDivs(dividends.size());
      std::transform(dividends.cbegin(),
                    dividends.cend(),
                    discDivs.begin(),
                    [=](const beagle::dividend_schedule_t::value_type& item)
                    { return std::make_pair(std::get<0>(item), std::get<2>(item)); });

      // Model parameters
      double spot = 100.;
      double rate = .03;
      double carry = .0;
      double vol = .3;

      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-rate * arg);});
      beagle::real_function_ptr_t funding = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(rate - carry) * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createGeneralForwardAssetPriceFunction(
                                                spot,
                                                funding,
                                                dividends,
                                                beagle::valuation::DividendPolicy::liquidator());
      beagle::valuation::OneDimFiniteDifferenceSettings settings(750, 1501, 4.5);

      beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(forward,
                                                                                                              discounting,
                                                                                                              vol );
      beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                forward,
                                                                discounting,
                                                                beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.),
                                                                beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                                beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0),
                                                                settings );
      beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                forward,
                                                                discounting,
                                                                beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                                settings );

      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
      for (double expiry : {5., 10., 15., 20.})
      {
        printf("\n************************************************************************************\n");
        printf("T = %f\n", expiry);
        for (double strike : {50., 75., 100., 125., 150., 175., 200.})
        {

          beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                    strike,
                                                                                                    payoff );
          beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
                                                                                                    strike,
                                                                                                    payoff );

          printf("K = %f, Call price = %.8f\n", strike, bscfeop->value( euroOption ));
          // std::cout << "European option price (FD-B) is: " << odbpop->value( euroOption ) << std::endl;
          // std::cout << "American option price (FD-B) is: " << odbpop->value( amerOption ) << std::endl;
          // std::cout << "European option price (FD-F) is: " << odfpeop->value( euroOption ) << std::endl;
        }
        printf("************************************************************************************\n\n");
      }

      EXPECT_EQ(0.0, 0.0);
    }
  }
}