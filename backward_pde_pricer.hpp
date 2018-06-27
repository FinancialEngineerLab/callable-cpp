#ifndef BACKWARD_PDE_PRICER_HPP
#define BACKWARD_PDE_PRICER_HPP

#include "pricer.hpp"
#include "dividend_policy.hpp"

namespace beagle
{
  namespace impl
  {
    struct OneDimensionalBackwardPDEOptionPricer : public Pricer
    {
      OneDimensionalBackwardPDEOptionPricer( double spot,
                                             double rate,
                                             const beagle::real_2d_function_ptr_t& volatility,
                                             const beagle::discrete_dividend_schedule_t& dividends,
                                             const beagle::dividend_policy_ptr_t& policy,
                                             int stepsPerAnnum );
      ~OneDimensionalBackwardPDEOptionPricer( void );
    public:
      virtual double optionValue( const beagle::option_ptr_t& option ) const override;
    private:
      using dbl_vec_t = std::vector<double>;
      using int_vec_t = std::vector<int>;
      void formLatticeForBackwardValuation( dbl_vec_t& times,
                                            dbl_vec_t&  ) const;
    private:
      double m_Spot;
      double m_Rate;
      beagle::real_2d_function_ptr_t m_Volatility;
      beagle::discrete_dividend_schedule_t m_Dividends;
      beagle::dividend_policy_ptr_t m_Policy;
      int m_StepsPerAnnum;
    };
  }
}

#endif