#ifndef BACKWARD_PDE_PRICER_HPP
#define BACKWARD_PDE_PRICER_HPP

#include "pricer.hpp"

namespace beagle
{
  namespace impl
  {
    struct BlackScholesBackwardPDEOptionPricer : public Pricer
    {
      BlackScholesBackwardPDEOptionPricer( double spot,
                                           double rate,
                                           double volatility,
                                           const beagle::discrete_dividend_schedule_t& dividends );
      ~BlackScholesBackwardPDEOptionPricer( void );
    public:
      virtual double optionValue( const beagle::option_ptr_t& option ) const override;
    private:
      double m_Spot;
      double m_Rate;
      double m_Volatility;
      beagle::discrete_dividend_schedule_t m_Dividends;
    };
  }
}

#endif