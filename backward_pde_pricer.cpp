#include "backward_pde_pricer.hpp"

namespace beagle
{

  namespace impl
  {
    BlackScholesBackwardPDEOptionPricer::BlackScholesBackwardPDEOptionPricer(
                                                   double spot,
                                                   double rate,
                                                   double volatility,
                                                   const beagle::discrete_dividend_schedule_t& dividends ) :
      m_Spot( spot ),
      m_Rate( rate ),
      m_Volatility( volatility ),
      m_Dividends( dividends )
    { }

    BlackScholesBackwardPDEOptionPricer::~BlackScholesBackwardPDEOptionPricer( void )
    { }

    double
    BlackScholesBackwardPDEOptionPricer::optionValue( const option_ptr_t& option ) const
    {
      return 0.0;
    }
  }

  beagle::pricer_ptr_t
  Pricer::formBlackScholesBackwardPDEOptionPricer( double spot,
                                                   double rate,
                                                   double volatility,
                                                   const beagle::discrete_dividend_schedule_t& dividends )
  {
    return std::make_shared<impl::BlackScholesBackwardPDEOptionPricer>( spot, rate, volatility, dividends );
  }
}