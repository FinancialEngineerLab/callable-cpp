#ifndef PRICER_HPP
#define PRICER_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace valuation
  {
    struct Pricer
    {
      Pricer( void );
      virtual ~Pricer( void );
    public:
      virtual double optionValue( const beagle::option_ptr_t& option ) const = 0;
    public:
      static beagle::pricer_ptr_t formBlackScholesClosedFormEuropeanOptionPricer(
                                                                          double spot,
                                                                          double rate,
                                                                          double volatility,
                                                                          const discrete_dividend_schedule_t& dividends );
      static beagle::pricer_ptr_t formOneDimensionalBackwardPDEOptionPricer(
                                                                    double spot,
                                                                     double rate,
                                                                     const beagle::real_2d_function_ptr_t& volatility,
                                                                     const beagle::discrete_dividend_schedule_t& dividends,
                                                                     const beagle::dividend_policy_ptr_t& policy,
                                                                     int stepsPerAnnum );
    };
  }
}


#endif