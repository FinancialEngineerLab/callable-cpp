#include <iostream>
#include "pricer.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"
#include "backward_pde_pricer.hpp"
#include "util.hpp"

int main( void )
{
  beagle::discrete_dividend_schedule_t dividends;
  // dividends.emplace_back(  0.5, 3. );
  // dividends.emplace_back(  1.5, 3. );
  // dividends.emplace_back(  2.5, 3. );
  // dividends.emplace_back(  3.5, 3. );
  // dividends.emplace_back(  4.5, 3. );
  // dividends.emplace_back(  5.5, 3. );
  // dividends.emplace_back(  6.5, 3. );
  // dividends.emplace_back(  7.5, 3. );
  // dividends.emplace_back(  8.5, 3. );
  // dividends.emplace_back(  9.5, 3. );
  // dividends.emplace_back( 10.5, 3. );
  // dividends.emplace_back( 11.5, 3. );
  // dividends.emplace_back( 12.5, 3. );
  // dividends.emplace_back( 13.5, 3. );
  // dividends.emplace_back( 14.5, 3. );
  // dividends.emplace_back( 15.5, 3. );
  // dividends.emplace_back( 16.5, 3. );
  // dividends.emplace_back( 17.5, 3. );
  // dividends.emplace_back( 18.5, 3. );
  // dividends.emplace_back( 19.5, 3. );
  // dividends.emplace_back( 20.5, 3. );
  // dividends.emplace_back( 21.5, 3. );
  // dividends.emplace_back( 22.5, 3. );
  // dividends.emplace_back( 23.5, 3. );

  try
  {
    beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer( 100., .03, .3, dividends );
    beagle::option_ptr_t euroOption = beagle::option::Option::createEuropeanOption( 5., 99.99, beagle::option::Payoff::put() );
    double value = bscfeop->optionValue( euroOption );
    std::cout << value << std::endl;
  }
  catch (const std::string& what)
  {
    std::cout << "A valuation error has occurred -- " << what << std::endl;
  }

  try
  {
    beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer( 100.,
                                                                                                         .03,
                                                                                                         beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(.3),
                                                                                                         1201,
                                                                                                         1601,
                                                                                                         5.5,
                                                                                                         dividends,
                                                                                                         beagle::valuation::DividendPolicy::liquidator(),
                                                                                                         beagle::math::InterpolationBuilder::linearWithFlatExtrapolation() );
    beagle::option_ptr_t option = beagle::option::Option::createAmericanOption( 5., 99.99, beagle::option::Payoff::put() );
    double value = bscfeop->optionValue( option );
    std::cout << value << std::endl;
  }
  catch (const std::string& what)
  {
    std::cout << "A valuation error has occurred -- " << what << std::endl;
  }

  return 0;
}