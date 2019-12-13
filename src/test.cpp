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
  std::cout << "\nStart of Test 1:\n";

  beagle::discrete_dividend_schedule_t dividends;
  dividends.emplace_back( 0.5, 6.0 );
  dividends.emplace_back( 1.5, 6.5 );
  dividends.emplace_back( 2.5, 7.0 );
  dividends.emplace_back( 3.5, 7.5 );
  dividends.emplace_back( 4.5, 8.0 );
  dividends.emplace_back( 5.5, 8.0 );
  dividends.emplace_back( 6.5, 8.0 );

  double expiry = 7.;
  double strike = 70.;
  beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();
  beagle::option_ptr_t euroOption = beagle::option::Option::createEuropeanOption( expiry,
                                                                                  strike,
                                                                                  payoff );
  beagle::option_ptr_t amerOption = beagle::option::Option::createAmericanOption( expiry,
                                                                                  strike,
                                                                                  payoff );

  beagle::valuation::FiniteDifferenceDetails fdDetails( 100., .06, .25, 1501, 1901, 7.5, dividends,
                                                        beagle::valuation::DividendPolicy::liquidator(),
                                                        beagle::math::InterpolationBuilder::linear());

  try
  {
    beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(fdDetails.spot(), 
                                                                                                             fdDetails.rate(), 
                                                                                                             fdDetails.volatility(), 
                                                                                                             fdDetails.dividends() );
    beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer(
                                                               fdDetails,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(fdDetails.volatility()) );
    beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                               fdDetails,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(fdDetails.volatility()) );

    std::cout << "European option price (CF)   is: " << bscfeop->optionValue( euroOption ) << std::endl;
    std::cout << "European option price (FD-B) is: " << odbpop->optionValue( euroOption ) << std::endl;
    std::cout << "American option price (FD-B) is: " << odbpop->optionValue( amerOption ) << std::endl;
    std::cout << "European option price (FD-F) is: " << odfpeop->optionValue( euroOption ) << std::endl;

    std::cout << "\nEnd of Test 1:\n";
  }
  catch (const std::string& what)
  {
    std::cout << "A valuation error has occurred -- " << what << std::endl;
  }
}

void test2( void )
{
  std::cout << "\nStart of Test 2:\n";

  beagle::discrete_dividend_schedule_t dividends;
  dividends.emplace_back( 0.4, 3.0 );
  dividends.emplace_back(0.6, 2.0);
  dividends.emplace_back(0.8, 1.0);

  beagle::valuation::FiniteDifferenceDetails fdDetails(100., .00, .25, 1001, 2001, 7.5, dividends,
                                                       beagle::valuation::DividendPolicy::liquidator(),
                                                       beagle::math::InterpolationBuilder::linear());

  beagle::dbl_vec_t expiries;
  beagle::dbl_vec_vec_t strikesColl;
  beagle::dbl_vec_vec_t pricesColl;
  beagle::test::generateEuropeanMarketQuotes(fdDetails, expiries, strikesColl, pricesColl );

  try
  {
    beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();
    beagle::pricer_ptr_t forwardPricer  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                 fdDetails,
                                                                 beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(fdDetails.volatility()) );

    beagle::real_2d_function_ptr_t localVol =
      beagle::math::RealTwoDimFunction::createBootstrappedLocalVolatilityFunction( expiries,
                                                                                   strikesColl,
                                                                                   pricesColl,
                                                                                   forwardPricer,
                                                                                   payoff );

    for (auto strike : strikesColl[0])
      std::cout << localVol->value( expiries[0], strike ) << std::endl;


    std::cout << "\nEnd of Test 2:\n";
  }
  catch (const std::string& what)
  {
    std::cout << what << std::endl;
  }
}

void test3( void )
{
  std::cout << "\nStart of Test 3:\n";

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
  beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();
  beagle::option_ptr_t euroOption = beagle::option::Option::createEuropeanOption( expiry,
                                                                                  strike,
                                                                                  payoff );
  beagle::option_ptr_t amerOption = beagle::option::Option::createAmericanOption( expiry,
                                                                                  strike,
                                                                                  payoff );

  double spot = 100.;
  double rate = .0;

  // CEV modeling
  double alpha = .3;
  double beta = .5;
  beagle::real_function_ptr_t localVolFunction = beagle::math::RealFunction::createUnaryFunction(
                                                              [alpha, beta](double x) { return alpha * std::pow(x, beta - 1.); } );
  beagle::real_2d_function_ptr_t localVolSurface =
            beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( beagle::dbl_vec_t(1U, 1.),
                                                                                    beagle::real_function_ptr_coll_t(1U, localVolFunction) );

  beagle::valuation::FiniteDifferenceDetails fdDetails(100., .00, .25, 1501, 1901, 7.5, dividends,
                                                       beagle::valuation::DividendPolicy::liquidator(),
                                                       beagle::math::InterpolationBuilder::linear());

  try
  {
    beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer(
                                                               fdDetails,
                                                               localVolSurface );
    beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                               fdDetails,
                                                               localVolSurface );

    std::cout << "European option price (FD-B) is: " << odbpop->optionValue( euroOption ) << std::endl;
    std::cout << "American option price (FD-B) is: " << odbpop->optionValue( amerOption ) << std::endl;
    std::cout << "European option price (FD-F) is: " << odfpeop->optionValue( euroOption ) << std::endl;

    std::cout << "\nEnd of Test 3:\n";
  }
  catch (const std::string& what)
  {
    std::cout << "A valuation error has occurred -- " << what << std::endl;
  }
}

int main( void )
{
  test1();
  test2();
  test3();

  return 0;
}
