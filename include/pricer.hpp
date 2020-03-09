#ifndef PRICER_HPP
#define PRICER_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace valuation
  {
    struct OneDimFiniteDifferenceSettings
    {
      OneDimFiniteDifferenceSettings( void );
      OneDimFiniteDifferenceSettings( int numTimeSteps,
                                      int numStateVariableSteps,
                                      double numGaussianStandardDeviations,
                                      const beagle::interp_builder_ptr_t& interp );
      OneDimFiniteDifferenceSettings( int numTimeSteps,
                                      int numStateVariableSteps,
                                      double numGaussianStandardDeviations,
                                      const beagle::interp_builder_ptr_t& interp,
                                      const beagle::dividend_policy_ptr_t& policy );
    //public:
    //  void formTimeSteps( double start,
    //                      double end,
    //                      beagle::dbl_vec_t& times,
    //                      beagle::int_vec_t& exDividendIndices ) const;
    //  void formStateVariableSteps(double expiry,
    //                              beagle::dbl_vec_t& logStateVariables,
    //                              beagle::dbl_vec_t& stateVariables) const;
    public:
      int numberOfTimeSteps( void ) const;
      int numberOfStateVariableSteps( void ) const;
      double numberOfGaussianStandardDeviations( void ) const;
      const beagle::dividend_policy_ptr_t& dividendPolicy( void ) const;
      const beagle::interp_builder_ptr_t& interpolationMethod( void ) const;
    private:
      int m_NumTimeSteps;
      int m_NumUnderlyingSteps;
      double m_NumStdev;
      beagle::dividend_policy_ptr_t m_Policy;
      beagle::interp_builder_ptr_t m_Interp;
    };

    struct FiniteDifferenceDetails
    {
      FiniteDifferenceDetails( void );
      FiniteDifferenceDetails(double spot,
                              double rate,
                              double volatility,
                              int stepsPerAnnum,
                              int stepsLogSpot,
                              double numStdev,
                              const beagle::discrete_dividend_schedule_t& dividends,
                              const beagle::dividend_policy_ptr_t& policy,
                              const beagle::interp_builder_ptr_t& interp);
    public:
      void formTimeSteps( double start,
                          double end,
                          beagle::dbl_vec_t& times,
                          beagle::int_vec_t& exDividendIndices ) const;
      void formStateVariableSteps(double expiry,
                                  beagle::dbl_vec_t& logStateVariables,
                                  beagle::dbl_vec_t& stateVariables) const;
    public:
      double spot() const;
      double rate() const;
      double volatility() const;
      int stepsPerAnnum() const;
      int numberOfStateVariableSteps() const;
      double numberOfStandardDeviations() const;
      const beagle::discrete_dividend_schedule_t& dividends() const;
      const beagle::dividend_policy_ptr_t& dividendPolicy() const;
      const beagle::interp_builder_ptr_t& interpolation() const;
    private:
      double m_Spot;
      double m_Rate;
      double m_Volatility;
      int m_StepsPerAnnum;
      int m_StepsLogSpot;
      double m_NumStdev;
      beagle::discrete_dividend_schedule_t m_Dividends;
      beagle::dividend_policy_ptr_t m_Policy;
      beagle::interp_builder_ptr_t m_Interp;
    };

    struct Pricer
    {
      virtual ~Pricer( void );
    public:
      virtual double value( const beagle::product_ptr_t& product ) const = 0;
    public:
      static beagle::pricer_ptr_t formBlackScholesClosedFormEuropeanOptionPricer(
                                                                          double spot,
                                                                          double rate,
                                                                          double volatility,
                                                                          const discrete_dividend_schedule_t& dividends );
      static beagle::pricer_ptr_t formOneDimensionalBackwardPDEOptionPricer(
                                                                     const FiniteDifferenceDetails& fdDetails,
                                                                     const beagle::real_2d_function_ptr_t& diffusion );
      static beagle::pricer_ptr_t formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                     const FiniteDifferenceDetails& fdDetails,
                                                                     const beagle::real_2d_function_ptr_t& volatility);
      static beagle::pricer_ptr_t formOneDimForwardPDEEuroOptionPricer(double spot,
                                                                       const beagle::real_2d_function_ptr_t& drift,
                                                                       const beagle::real_2d_function_ptr_t& volatility,
                                                                       const beagle::real_function_ptr_t& rate,
                                                                       const beagle::discrete_dividend_schedule_t& dividends,
                                                                       const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
      static beagle::pricer_ptr_t formOneDimBackwardPDEOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                    const beagle::real_function_ptr_t& discounting,
                                                                    const beagle::real_2d_function_ptr_t& drift,
                                                                    const beagle::real_2d_function_ptr_t& volatility,
                                                                    const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
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
        virtual const FiniteDifferenceDetails& finiteDifferenceDetails(void) const = 0;
      };

      struct CloneWithNewLocalVolatilitySurface
      {
        virtual ~CloneWithNewLocalVolatilitySurface( void );
      public:
        virtual beagle::pricer_ptr_t createPricerWithNewLocalVolatilitySurface( const beagle::real_2d_function_ptr_t& vol ) const = 0;
      };

      struct CloneWithNewModelParameters
      {
        virtual ~CloneWithNewModelParameters( void );
      public:
        virtual beagle::pricer_ptr_t createPricerWithNewModelParameters( const beagle::dbl_vec_t& parameters ) const = 0;
      };
    }
  }
}


#endif
