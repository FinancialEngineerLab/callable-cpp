#ifndef EURO_PRICER_HPP
#define EURO_PRICER_HPP

#include "pricer.hpp"

namespace beagle
{
  namespace impl
  {  
    struct BlackScholesClosedFormEuropeanOptionPricer : public Pricer
    {
      BlackScholesClosedFormEuropeanOptionPricer( double spot,
                                                  double rate,
                                                  double volatility,
                                                  const beagle::discrete_dividend_schedule_t& dividends );
      ~BlackScholesClosedFormEuropeanOptionPricer( void );
    public:
      virtual double optionValue( const beagle::option_ptr_t& option ) const override;
    private:
      void calculateAdjustedSpotAndStrike( double expiry,
                                           double strike ) const;
    private:
      double m_Spot;
      double m_Rate;
      double m_Volatility;
      beagle::discrete_dividend_schedule_t m_Dividends;

      mutable double m_AdjustedSpot;
      mutable double m_AdjustedStrike;
    };
  }
}



#endif