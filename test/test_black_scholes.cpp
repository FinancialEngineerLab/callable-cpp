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
    TEST(test_black_scholes, AmericanPricing)
    {
      // Model parameters
      double spot = 45.;
      double rate = .02;
      double carry = .01;
      double vol = .2;

      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-rate * arg);});
      beagle::real_function_ptr_t funding = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(rate - carry) * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                funding);
      beagle::valuation::OneDimFiniteDifferenceSettings settings(750, 1501, 4.5);
      beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                forward,
                                                                discounting,
                                                                beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.),
                                                                beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                                beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0),
                                                                settings);

      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
      double expiry = 1.;
      double strike = 45.;
      beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
                                                                                                strike,
                                                                                                payoff );
      double value = odbpop->value( amerOption );
      std::cout << value << std::endl;
    }

    TEST(test_black_scholes, AnnualDividends)
    {
      // Set up dividend
      beagle::dividend_schedule_t dividends;
      dividends.reserve(20);
      for (int i = 0; i < 20; ++i)
        dividends.emplace_back(i + .5, 0.0, 3.0);

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
      // beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
      //                                                           forward,
      //                                                           discounting,
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.),
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0),
      //                                                           settings );
      // beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
      //                                                           forward,
      //                                                           discounting,
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
      //                                                           settings );

      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();

      // Test results are taken from https://arxiv.org/pdf/1008.3880.pdf
      beagle::dbl_vec_t expiries = {5., 10., 15., 20.};
      beagle::dbl_vec_t strikes = {50., 75., 100., 125., 150., 175., 200.};
      beagle::dbl_mat_t prices = {{47.14, 33.85, 24.42, 17.79, 13.12,  9.79,  7.39},
                                  {46.85, 38.21, 31.66, 26.58, 22.56, 19.34, 16.71},
                                  {46.49, 40.49, 35.73, 31.85, 28.63, 25.91, 23.59},
                                  {46.10, 41.76, 38.23, 35.26, 32.71, 30.50, 28.56}};
      for (int i = 0; i < expiries.size(); ++i)
      {
        double expiry = expiries[i];

        // printf("\n************************************************************************************\n");
        // printf("T = %3.0f\n", expiry);
        for (int j = 0; j < strikes.size(); ++j)
        {
          double strike = strikes[j];
          beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                    strike,
                                                                                                    payoff );
          double value = bscfeop->value( euroOption );
          EXPECT_NEAR(value / prices[i][j], 1., 5e-4);

          // printf("K = %3.0f, Call price = %.8f\n", strike, value);
          // beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
          //                                                                                           strike,
          //                                                                                           payoff );
          // std::cout << "European option price (FD-B) is: " << odbpop->value( euroOption ) << std::endl;
          // std::cout << "American option price (FD-B) is: " << odbpop->value( amerOption ) << std::endl;
          // std::cout << "European option price (FD-F) is: " << odfpeop->value( euroOption ) << std::endl;
        }
        // printf("************************************************************************************\n\n");
      }
    }

    TEST(test_black_scholes, WeeklyDividends)
    {
      // Set up dividend
      beagle::dividend_schedule_t dividends;
      dividends.reserve(1043);
      for (int i = 0; i < 1043; ++i)
        dividends.emplace_back((i + 1) * 7. / 365., 0.0, 2.0);

      beagle::discrete_dividend_schedule_t discDivs(dividends.size());
      std::transform(dividends.cbegin(),
                     dividends.cend(),
                     discDivs.begin(),
                     [=](const beagle::dividend_schedule_t::value_type& item)
                     { return std::make_pair(std::get<0>(item), std::get<2>(item)); });

      // Model parameters
      double spot = 3000.;
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
      // beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
      //                                                           forward,
      //                                                           discounting,
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.),
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0),
      //                                                           settings );
      // beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
      //                                                           forward,
      //                                                           discounting,
      //                                                           beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
      //                                                           settings );

      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();

      // Test results are taken from https://arxiv.org/pdf/1008.3880.pdf
      beagle::dbl_vec_t expiries = {5., 10., 15., 20.};
      beagle::dbl_vec_t strikes = {1500., 2250., 3000., 3750., 4500., 5250., 6000.};
      beagle::dbl_mat_t prices = {{1359.87,  972.69,  699.68, 508.73, 374.47, 279.07, 210.47},
                                  {1319.68, 1075.04,  889.96, 746.72, 633.69, 543.02, 469.25},
                                  {1288.47, 1122.33,  990.66, 883.42, 794.31, 719.10, 654.79},
                                  {1264.53, 1145.94, 1049.44, 968.43, 899.04, 838.71, 785.66}};
      for (int i = 0; i < expiries.size(); ++i)
      {
        double expiry = expiries[i];

        // printf("\n************************************************************************************\n");
        // printf("T = %3.0f\n", expiry);
        for (int j = 0; j < strikes.size(); ++j)
        {
          double strike = strikes[j];
          beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                    strike,
                                                                                                    payoff );
          double value = bscfeop->value( euroOption );
          EXPECT_NEAR(value / prices[i][j], 1., 5e-3);

          // printf("K = %3.0f, Call price = %.8f\n", strike, value);
          // beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
          //                                                                                           strike,
          //                                                                                           payoff );
          // std::cout << "European option price (FD-B) is: " << odbpop->value( euroOption ) << std::endl;
          // std::cout << "American option price (FD-B) is: " << odbpop->value( amerOption ) << std::endl;
          // std::cout << "European option price (FD-F) is: " << odfpeop->value( euroOption ) << std::endl;
        }
        // printf("************************************************************************************\n\n");
      }
    }
  }
}
