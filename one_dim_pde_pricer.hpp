#ifndef ONE_DIM_PDE_PRICER_HPP
#define ONE_DIM_PDE_PRICER_HPP

#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "util.hpp"
#include "dividend_policy.hpp"
#include "real_function.hpp"
#include "interpolation_builder.hpp"

// #include <fstream>
// std::ofstream out("interpolation.txt");

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct OneDimensionalPDEOptionPricer : public Pricer,
                                             public beagle::valuation::mixins::FiniteDifference
      {
        using two_dbl_t = std::pair<double, double>;

        OneDimensionalPDEOptionPricer( double spot,
                                       double rate,
                                       const beagle::real_2d_function_ptr_t& volatility,
                                       int stepsPerAnnum,
                                       int stepsLogSpot,
                                       double numStdev,
                                       const beagle::discrete_dividend_schedule_t& dividends,
                                       const beagle::dividend_policy_ptr_t& policy,
                                       const beagle::interp_builder_ptr_t& interp ) :
          m_Spot( spot ),
          m_Rate( rate ),
          m_Volatility( volatility ),
          m_StepsPerAnnum( stepsPerAnnum ),
          m_StepsLogSpot( stepsLogSpot ),
          m_NumStdev( numStdev ),
          m_Dividends( dividends ),
          m_Policy( policy ),
          m_Interp( interp )
        { }
        ~OneDimensionalPDEOptionPricer( void )
        { }
      public:
        virtual double spot( void ) const
        {
          return m_Spot;
        }
        virtual double rate( void ) const
        {
          return m_Rate;
        }
        virtual const beagle::real_2d_function_ptr_t& volatility( void ) const
        {
          return m_Volatility;
        }
        virtual int stepsPerAnnum( void ) const
        {
          return m_StepsPerAnnum;
        }
        virtual int numberOfStateVariableSteps( void ) const
        {
          return m_StepsLogSpot;
        }
        virtual int numberOfStandardDeviations( void ) const
        {
          return m_NumStdev;
        }
        virtual const beagle::discrete_dividend_schedule_t& dividends( void ) const
        {
          return m_Dividends;
        }
        virtual const beagle::dividend_policy_ptr_t& dividendPolicy( void ) const
        {
          return m_Policy;
        }
        virtual const interp_builder_ptr_t& interpolation( void ) const
        {
          return m_Interp;
        }
      private:
        double m_Spot;
        double m_Rate;
        beagle::real_2d_function_ptr_t m_Volatility;
        int m_StepsPerAnnum;
        int m_StepsLogSpot;
        double m_NumStdev;
        beagle::discrete_dividend_schedule_t m_Dividends;
        beagle::dividend_policy_ptr_t m_Policy;
        beagle::interp_builder_ptr_t m_Interp;
      };
    }
  }
}

#endif