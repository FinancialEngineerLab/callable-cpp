#include <iostream>
#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"
#include "util.hpp"
#include "round_trip.hpp"

void test1( void )
{
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

  double spot = 100.;
  double rate = .06;
  double vol  = .25;

  try
  {
    beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer( spot, rate, vol, dividends );
    beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer(
                                                               spot,
                                                               rate,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                               1501,
                                                               1901,
                                                               7.5,
                                                               dividends,
                                                               beagle::valuation::DividendPolicy::liquidator(),
                                                               beagle::math::InterpolationBuilder::linear() );
    beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                               spot,
                                                               rate,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                               1501,
                                                               1901,
                                                               7.5,
                                                               dividends,
                                                               beagle::valuation::DividendPolicy::liquidator(),
                                                               beagle::math::InterpolationBuilder::linear() );
    std::cout << "European option price (CF)   is: " << bscfeop->optionValue( euroOption ) << std::endl;
    std::cout << "European option price (FD-B) is: " << odbpop->optionValue( euroOption ) << std::endl;
    std::cout << "American option price (FD-B) is: " << odbpop->optionValue( amerOption ) << std::endl;
    std::cout << "European option price (FD-F) is: " << odfpeop->optionValue( euroOption ) << std::endl;
  }
  catch (const std::string& what)
  {
    std::cout << "A valuation error has occurred -- " << what << std::endl;
  }
}

void test2( void )
{
  beagle::dbl_vec_t expiries;
  beagle::dbl_vec_vec_t strikesColl;
  beagle::dbl_vec_vec_t pricesColl;
  beagle::test::generateEuropeanMarketQuotes( expiries, strikesColl, pricesColl );

  double spot = 100.;
  double rate = .00;
  double vol  = .25;
  try
  {
    beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();
    beagle::pricer_ptr_t forwardPricer  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                 spot,
                                                                 rate,
                                                                 beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                                 1001,
                                                                 1001,
                                                                 5.,
                                                                 beagle::discrete_dividend_schedule_t(),
                                                                 beagle::valuation::DividendPolicy::liquidator(),
                                                                 beagle::math::InterpolationBuilder::linear() );

    beagle::real_2d_function_ptr_t localVol =
      beagle::math::RealTwoDimFunction::createBootstrappedLocalVolatilityFunction( expiries,
                                                                                   strikesColl,
                                                                                   pricesColl,
                                                                                   forwardPricer,
                                                                                   payoff );

    for (auto strike : strikesColl[0])
      std::cout << localVol->value( expiries[0], strike ) << std::endl;
  }
  catch (const std::string& what)
  {
    std::cout << what << std::endl;
  }
}

void test3( void )
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
  double strike = 90.;
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

  try
  {
    beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer(
                                                               spot,
                                                               rate,
                                                               localVolSurface,
                                                               1501,
                                                               1901,
                                                               7.5,
                                                               dividends,
                                                               beagle::valuation::DividendPolicy::liquidator(),
                                                               beagle::math::InterpolationBuilder::linear() );
    beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                               spot,
                                                               rate,
                                                               localVolSurface,
                                                               1501,
                                                               2901,
                                                               9.5,
                                                               dividends,
                                                               beagle::valuation::DividendPolicy::liquidator(),
                                                               beagle::math::InterpolationBuilder::linear() );
    std::cout << "European option price (FD-B) is: " << odbpop->optionValue( euroOption ) << std::endl;
    std::cout << "American option price (FD-B) is: " << odbpop->optionValue( amerOption ) << std::endl;
    std::cout << "European option price (FD-F) is: " << odfpeop->optionValue( euroOption ) << std::endl;
  }
  catch (const std::string& what)
  {
    std::cout << "A valuation error has occurred -- " << what << std::endl;
  }
}

int main( void )
{
  test2();

  return 0;
}
