#include <iostream>
#include "euro_pricer.hpp"

int main( void )
{
  beagle::discrete_dividend_schedule_t dividends;
  dividends.emplace_back(  0.5, 3. );
  dividends.emplace_back(  1.5, 3. );
  dividends.emplace_back(  2.5, 3. );
  dividends.emplace_back(  3.5, 3. );
  dividends.emplace_back(  4.5, 3. );
  dividends.emplace_back(  5.5, 3. );
  dividends.emplace_back(  6.5, 3. );
  dividends.emplace_back(  7.5, 3. );
  dividends.emplace_back(  8.5, 3. );
  dividends.emplace_back(  9.5, 3. );
  dividends.emplace_back( 10.5, 3. );
  dividends.emplace_back( 11.5, 3. );
  dividends.emplace_back( 12.5, 3. );
  dividends.emplace_back( 13.5, 3. );
  dividends.emplace_back( 14.5, 3. );
  dividends.emplace_back( 15.5, 3. );
  dividends.emplace_back( 16.5, 3. );
  dividends.emplace_back( 17.5, 3. );
  dividends.emplace_back( 18.5, 3. );
  dividends.emplace_back( 19.5, 3. );
  dividends.emplace_back( 20.5, 3. );
  dividends.emplace_back( 21.5, 3. );
  dividends.emplace_back( 22.5, 3. );
  dividends.emplace_back( 23.5, 3. );

  beagle::EuropeanOptionClosedFormPricer opt( 100,
                                              .03,
                                              .3,
                                              dividends );
  double value = opt.callOptionValue( 20, 100 );

  std::cout << value << std::endl;
  return 0;
}