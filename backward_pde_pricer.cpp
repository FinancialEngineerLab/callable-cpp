#include "backward_pde_pricer.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      OneDimensionalBackwardPDEOptionPricer::OneDimensionalBackwardPDEOptionPricer(
                                                     double spot,
                                                     double rate,
                                                     const beagle::real_2d_function_ptr_t& volatility,
                                                     const beagle::discrete_dividend_schedule_t& dividends,
                                                     const beagle::dividend_policy_ptr_t& policy,
                                                     int stepsPerAnnum ) :
        m_Spot( spot ),
        m_Rate( rate ),
        m_Volatility( volatility ),
        m_Dividends( dividends ),
        m_Policy( policy ),
        m_StepsPerAnnum( stepsPerAnnum )
      { }

      OneDimensionalBackwardPDEOptionPricer::~OneDimensionalBackwardPDEOptionPricer( void )
      { }

      double
      OneDimensionalBackwardPDEOptionPricer::optionValue( const beagle::option_ptr_t& option ) const
      {
        return 0.0;
      }
    }

    beagle::pricer_ptr_t
    Pricer::formOneDimensionalBackwardPDEOptionPricer( double spot,
                                                       double rate,
                                                       const beagle::real_2d_function_ptr_t& volatility,
                                                       const beagle::discrete_dividend_schedule_t& dividends,
                                                       const beagle::dividend_policy_ptr_t& policy,
                                                       int stepsPerAnnum )
    {
      return std::make_shared<impl::OneDimensionalBackwardPDEOptionPricer>( spot, rate, volatility, dividends, policy, stepsPerAnnum );
    }
  }
}