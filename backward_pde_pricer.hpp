#ifndef BACKWARD_PDE_PRICER_HPP
#define BACKWARD_PDE_PRICER_HPP

#include "pricer.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      using two_dbl_t = std::pair<double, double>;
      struct OneDimensionalBackwardPDEOptionPricer : public Pricer
      {
        OneDimensionalBackwardPDEOptionPricer( double spot,
                                               double rate,
                                               const beagle::real_2d_function_ptr_t& volatility,
                                               int stepsPerAnnum,
                                               int stepsLogSpot,
                                               double numStdev,
                                               const beagle::discrete_dividend_schedule_t& dividends,
                                               const beagle::dividend_policy_ptr_t& policy,
                                               const beagle::interp_builder_ptr_t& interp );
        ~OneDimensionalBackwardPDEOptionPricer( void );
      public:
        virtual double optionValue( const beagle::option_ptr_t& option ) const override;
      private:
        void formLatticeForBackwardValuation( double expiry,
                                              beagle::dbl_vec_t& times,
                                              beagle::dbl_vec_t& logSpot,
                                              beagle::int_vec_t& exDividendIndices ) const;
        two_dbl_t boundaryCondition( const beagle::payoff_ptr_t& payoff,
                                     const two_dbl_t& boundarySpots,
                                     double strike,
                                     double time,
                                     bool isAmerican ) const;
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