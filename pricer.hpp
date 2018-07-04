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
                                                                     int stepsPerAnnum,
                                                                     int stepsLogSpot,
                                                                     double numStdev,
                                                                     const beagle::discrete_dividend_schedule_t& dividends,
                                                                     const beagle::dividend_policy_ptr_t& policy,
                                                                     const beagle::interp_builder_ptr_t& interp );
      static beagle::pricer_ptr_t formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                     double spot,
                                                                     double rate,
                                                                     const beagle::real_2d_function_ptr_t& volatility,
                                                                     int stepsPerAnnum,
                                                                     int stepsLogSpot,
                                                                     double numStdev,
                                                                     const beagle::discrete_dividend_schedule_t& dividends,
                                                                     const beagle::dividend_policy_ptr_t& policy,
                                                                     const beagle::interp_builder_ptr_t& interp );
    };

    namespace mixins
    {
      struct OptionValueCollectionProvider
      {
        virtual ~OptionValueCollectionProvider( void );
      public:
        virtual void formInitialOptionValueCollection( const beagle::payoff_ptr_t& payoff,
                                                       const beagle::dbl_vec_t& strikes,
                                                       beagle::dbl_vec_t& prices ) const = 0;
        virtual void optionValueCollection( double start,
                                            double end,
                                            const beagle::payoff_ptr_t& payoff,
                                            const beagle::dbl_vec_t& logStrikes,
                                            const beagle::dbl_vec_t& strikes,
                                            beagle::dbl_vec_t& prices ) const = 0;
      };

      struct FiniteDifference
      {
        virtual ~FiniteDifference( void );
      public:
        virtual void formTimeSteps( double start,
                                    double end,
                                    beagle::dbl_vec_t& times,
                                    beagle::int_vec_t& exDividendIndices ) const;
        virtual void formStateVariableSteps( double expiry,
                                             beagle::dbl_vec_t& logStateVariables,
                                             beagle::dbl_vec_t& stateVariables ) const;
      protected:
        virtual int stepsPerAnnum( void ) const = 0;
        virtual const beagle::discrete_dividend_schedule_t& dividends( void ) const = 0;
        virtual double spot( void ) const = 0;
        virtual double rate( void ) const = 0;
        virtual const beagle::real_2d_function_ptr_t& volatility( void ) const = 0;
        virtual int numberOfStandardDeviations( void ) const = 0;
        virtual int numberOfStateVariableSteps( void ) const = 0;
      };
    }
  }
}


#endif