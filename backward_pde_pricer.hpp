#ifndef BACKWARD_PDE_PRICER_HPP
#define BACKWARD_PDE_PRICER_HPP

#include "pricer.hpp"

namespace beagle
{
  namespace impl
  {
    struct OneDimensionalBackwardPDEOptionPricer : public Pricer
    {
      OneDimensionalBackwardPDEOptionPricer( double spot,
                                             double rate,
                                             const real_2d_function_ptr_t& volatility,
                                             const beagle::discrete_dividend_schedule_t& dividends );
      ~OneDimensionalBackwardPDEOptionPricer( void );
    public:
      virtual double optionValue( const beagle::option_ptr_t& option ) const override;
    private:
      double m_Spot;
      double m_Rate;
      real_2d_function_ptr_t m_Volatility;
      beagle::discrete_dividend_schedule_t m_Dividends;
    };
  }
}

#endif