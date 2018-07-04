#include <iostream>
#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"
#include "util.hpp"

int main( void )
{
  beagle::discrete_dividend_schedule_t dividends;
  dividends.emplace_back( 0.5, 6.0 );
  dividends.emplace_back( 1.5, 6.5 );
  // dividends.emplace_back( 2.5, 7.0 );
  // dividends.emplace_back( 3.5, 7.5 );
  // dividends.emplace_back( 4.5, 8.0 );
  // dividends.emplace_back( 5.5, 8.0 );
  // dividends.emplace_back( 6.5, 8.0 );

  double expiry = 2.;
  double strike = 100.;
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

  return 0;
}