#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"
#include "util.hpp"
#include "round_trip.hpp"

#include <iostream>

void test1( void )
{
  std::cout << "\nStart of Test 1:\n\n";

  beagle::discrete_dividend_schedule_t dividends;
  //dividends.emplace_back( 0.5, 3.0 );
  //dividends.emplace_back( 1.5, 3.0 );
  //dividends.emplace_back( 2.5, 3.0 );
  //dividends.emplace_back( 3.5, 3.0 );
  //dividends.emplace_back( 4.5, 3.0 );
  //dividends.emplace_back( 5.5, 8.0 );
  //dividends.emplace_back( 6.5, 8.0 );

  double expiry = .05;
  double strike = 100.;
  beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::put();
  beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                            strike,
                                                                                            payoff );
  beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
                                                                                            strike,
                                                                                            payoff );

  beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                            [](double arg) { return std::exp(-.3 * arg);});
  beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                            100, discounting);

  beagle::valuation::FiniteDifferenceDetails fdDetails( 100., .3, .3, 1501, 1901, 7.5, dividends,
                                                        beagle::valuation::DividendPolicy::liquidator(),
                                                        beagle::math::InterpolationBuilder::linear());

  try
  {
    beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(fdDetails.spot(), 
                                                                                                             fdDetails.rate(), 
                                                                                                             fdDetails.volatility(), 
                                                                                                             fdDetails.dividends() );
    beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                               forward,
                                                               discounting,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.),
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.3),
                                                               beagle::valuation::OneDimFiniteDifferenceSettings(1501, 1901, 7.5) );
    beagle::pricer_ptr_t odbpop2  = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer(
                                                               fdDetails,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(fdDetails.volatility()) );
    beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                               fdDetails,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(fdDetails.volatility()) );

    std::cout << "European option price (CF)   is: " << bscfeop->value( euroOption ) << std::endl;
    std::cout << "European option price (FD-B) is: " << odbpop->value( euroOption ) << std::endl;
    std::cout << "European option price (FD-B) is: " << odbpop2->value( euroOption ) << std::endl;
    std::cout << "American option price (FD-B) is: " << odbpop->value( amerOption ) << std::endl;
    std::cout << "American option price (FD-B) is: " << odbpop2->value( amerOption ) << std::endl;
    //std::cout << "European option price (FD-F) is: " << odfpeop->value( euroOption ) << std::endl;

    std::cout << "\nEnd of Test 1\n";
  }
  catch (const std::string& what)
  {
    std::cout << "A valuation error has occurred -- " << what << std::endl;
  }
}

void test2( void )
{
  std::cout << "\nStart of Test 2:\n\n";

  beagle::discrete_dividend_schedule_t dividends;
  // dividends.emplace_back(0.4, 3.0);
  // dividends.emplace_back(0.6, 2.0);
  // dividends.emplace_back(0.8, 1.0);

  beagle::valuation::FiniteDifferenceDetails fdDetails(100., .00, .5, 500, 1500, 5.5, dividends,
                                                       beagle::valuation::DividendPolicy::liquidator(),
                                                       beagle::math::InterpolationBuilder::linear());

  beagle::dbl_vec_t expiries;
  beagle::dbl_vec_vec_t strikesColl;
  beagle::dbl_vec_vec_t pricesColl;
  beagle::test::generateEuropeanMarketQuotes(fdDetails, expiries, strikesColl, pricesColl );

  beagle::dbl_vec_t initialGuesses(expiries.size(), .6);

  try
  {
    beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
    beagle::interp_builder_ptr_t interp = beagle::math::InterpolationBuilder::linear();

    beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                 fdDetails,
                                                                 beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(fdDetails.volatility()) );

    beagle::real_2d_function_ptr_t localVol =
      beagle::math::RealTwoDimFunction::createBootstrappedLocalVolatilityFunction( expiries,
                                                                                   initialGuesses,
                                                                                   strikesColl,
                                                                                   pricesColl,
                                                                                   forwardPricer,
                                                                                   payoff,
                                                                                   interp );

    // Output local volatility surface
    for (beagle::dbl_vec_t::size_type i = 0; i < expiries.size() - 1; ++i)
    {
      for (beagle::dbl_vec_t::size_type j = 0; j < strikesColl[i].size(); ++j)
        std::cout << expiries[i] - .1 << "\t"
                  << strikesColl[i][j] << "\t"
                  << localVol->value(expiries[i] - .1, strikesColl[i][j]) << "\n";

      std::cout << "\n";
    }

    std::cout << "\n";

    // Output price quotes
    for (beagle::dbl_vec_t::size_type i = 0; i < expiries.size() - 1; ++i)
    {
      for (beagle::dbl_vec_t::size_type j = 0; j < strikesColl[i].size(); ++j)
        std::cout << pricesColl[i][j] << "\t";

      std::cout << "\n";
    }

    std::cout << "\n";

    // Output calibrated prices
    beagle::pricer_ptr_t calibratedPricer = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer( fdDetails, localVol );
    for (beagle::dbl_vec_t::size_type i = 0; i < expiries.size() - 1; ++i)
    {
      for (beagle::dbl_vec_t::size_type j = 0; j < strikesColl[i].size(); ++j)
      {
        beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption(expiries[i],
                                                                                                 strikesColl[i][j],
                                                                                                 payoff);
        std::cout << calibratedPricer->value(euroOption) << "\t";
      }

      std::cout << "\n";
    }

    std::cout << "\n";

    std::cout << "\nEnd of Test 2\n";
  }
  catch (const std::string& what)
  {
    std::cout << what << std::endl;
  }
}
//
//void test3( void )
//{
//  std::cout << "\nStart of Test 3:\n\n";
//
//  beagle::discrete_dividend_schedule_t dividends;
//  // dividends.emplace_back( 0.5, 6.0 );
//  // dividends.emplace_back( 1.5, 6.5 );
//  // dividends.emplace_back( 2.5, 7.0 );
//  // dividends.emplace_back( 3.5, 7.5 );
//  // dividends.emplace_back( 4.5, 8.0 );
//  // dividends.emplace_back( 5.5, 8.0 );
//  // dividends.emplace_back( 6.5, 8.0 );
//
//  double expiry = 1.;
//  double strike = 100.;
//  beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();
//  beagle::option_ptr_t euroOption = beagle::option::Option::createEuropeanOption( expiry,
//                                                                                  strike,
//                                                                                  payoff );
//  beagle::option_ptr_t amerOption = beagle::option::Option::createAmericanOption( expiry,
//                                                                                  strike,
//                                                                                  payoff );
//
//  double spot = 100.;
//  double rate = .0;
//
//  // CEV modeling
//  double alpha = .3;
//  double beta = .5;
//  beagle::real_function_ptr_t localVolFunction = beagle::math::RealFunction::createUnaryFunction(
//                                                              [alpha, beta](double x) { return alpha * std::pow(x, beta - 1.); } );
//  beagle::real_2d_function_ptr_t localVolSurface =
//            beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( beagle::dbl_vec_t(1U, 1.),
//                                                                                    beagle::real_function_ptr_coll_t(1U, localVolFunction) );
//
//  beagle::valuation::FiniteDifferenceDetails fdDetails(100., .00, .25, 1501, 1901, 7.5, dividends,
//                                                       beagle::valuation::DividendPolicy::liquidator(),
//                                                       beagle::math::InterpolationBuilder::linear());
//
//  try
//  {
//    beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer(
//                                                               fdDetails,
//                                                               localVolSurface );
//    beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
//                                                               fdDetails,
//                                                               localVolSurface );
//
//    std::cout << "European option price (FD-B) is: " << odbpop->optionValue( euroOption ) << std::endl;
//    std::cout << "American option price (FD-B) is: " << odbpop->optionValue( amerOption ) << std::endl;
//    std::cout << "European option price (FD-F) is: " << odfpeop->optionValue( euroOption ) << std::endl;
//
//    std::cout << "\nEnd of Test 3\n";
//  }
//  catch (const std::string& what)
//  {
//    std::cout << "A valuation error has occurred -- " << what << std::endl;
//  }
//}

void test4(void)
{
  std::cout << "\nStart of Test 4:\n\n";

  beagle::dbl_vec_t diag{ -2.6, -2.6, -2.6, -2.6 };
  beagle::dbl_vec_t upper{ 1., 1., 1., 1. };
  beagle::dbl_vec_t lower{ 1., 1., 1., 1. };
  beagle::dbl_vec_t rhs{ -240., 0., 0., -150. };

  beagle::util::tridiagonalSolve(rhs, diag, upper, lower);
  for (auto result : rhs)
    std::cout << result << "\n";

  std::cout << "\nEnd of Test 4:\n\n";
}

int main( void )
{
  test1();
  //test2();
  //test3();
  //test4();

  return 0;
}
