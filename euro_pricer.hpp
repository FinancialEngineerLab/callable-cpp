#ifndef EURO_PRICER_HPP
#define EURO_PRICER_HPP

#include <vector>

namespace beagle
{
  using discrete_dividend_schedule_t = std::vector< std::pair<double, double> >;

  struct EuropeanOptionClosedFormPricer
  {
    EuropeanOptionClosedFormPricer( double spot,
                                    double rate,
                                    double volatility,
                                    const discrete_dividend_schedule_t& dividends );
    ~EuropeanOptionClosedFormPricer( void );
  public:
    double callOptionValue( double expiry,
                            double strike ) const;
  private:
    double m_Spot;
    double m_Rate;
    double m_Volatility;
    discrete_dividend_schedule_t m_Dividends;
  };

}



#endif