#ifndef PRICER_HPP
#define PRICER_HPP

#include <memory>
#include <vector>
#include "option.hpp"
#include "real_2d_function.hpp"

namespace beagle
{
  struct Pricer;
  using pricer_ptr_t = std::shared_ptr<Pricer>;
  using discrete_dividend_schedule_t = std::vector< std::pair<double, double> >;

  struct Pricer
  {
    Pricer( void );
    virtual ~Pricer( void );
  public:
    virtual double optionValue( const beagle::option_ptr_t& option ) const = 0;
  public:
    static pricer_ptr_t formBlackScholesClosedFormEuropeanOptionPricer( double spot,
                                                                        double rate,
                                                                        double volatility,
                                                                        const discrete_dividend_schedule_t& dividends );
    static pricer_ptr_t formOneDimensionalBackwardPDEOptionPricer( double spot,
                                                                   double rate,
                                                                   const real_2d_function_ptr_t& volatility,
                                                                   const discrete_dividend_schedule_t& dividends );
  };
}


#endif