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
                                      beagle::dbl_t numGaussianStandardDeviations );
      OneDimFiniteDifferenceSettings( int numTimeSteps,
                                      int numStateVariableSteps,
                                      beagle::dbl_t numGaussianStandardDeviations,
                                      const beagle::interp_builder_ptr_t& interp,
                                      const beagle::dividend_policy_ptr_t& policy );
    //public:
    //  void formTimeSteps( beagle::dbl_t start,
    //                      beagle::dbl_t end,
    //                      beagle::dbl_vec_t& times,
    //                      beagle::int_vec_t& exDividendIndices ) const;
    //  void formStateVariableSteps(beagle::dbl_t expiry,
    //                              beagle::dbl_vec_t& logStateVariables,
    //                              beagle::dbl_vec_t& stateVariables) const;
    public:
      int numberOfTimeSteps( void ) const;
      int numberOfStateVariableSteps( void ) const;
      beagle::dbl_t numberOfGaussianStandardDeviations( void ) const;
      const beagle::dividend_policy_ptr_t& dividendPolicy( void ) const;
      const beagle::interp_builder_ptr_t& interpolationMethod( void ) const;
    private:
      int m_NumTimeSteps;
      int m_NumUnderlyingSteps;
      beagle::dbl_t m_NumStdev;
      beagle::dividend_policy_ptr_t m_Policy;
      beagle::interp_builder_ptr_t m_Interp;
    };

    struct FiniteDifferenceDetails
    {
      FiniteDifferenceDetails( void );
      FiniteDifferenceDetails(beagle::dbl_t spot,
                              beagle::dbl_t rate,
                              beagle::dbl_t volatility,
                              int stepsPerAnnum,
                              int stepsLogSpot,
                              beagle::dbl_t numStdev,
                              const beagle::discrete_dividend_schedule_t& dividends,
                              const beagle::dividend_policy_ptr_t& policy,
                              const beagle::interp_builder_ptr_t& interp);
    public:
      void formTimeSteps( beagle::dbl_t start,
                          beagle::dbl_t end,
                          beagle::dbl_vec_t& times,
                          beagle::int_vec_t& exDividendIndices ) const;
      void formStateVariableSteps(beagle::dbl_t expiry,
                                  beagle::dbl_vec_t& logStateVariables,
                                  beagle::dbl_vec_t& stateVariables) const;
    public:
      beagle::dbl_t spot() const;
      beagle::dbl_t rate() const;
      beagle::dbl_t volatility() const;
      int stepsPerAnnum() const;
      int numberOfStateVariableSteps() const;
      beagle::dbl_t numberOfStandardDeviations() const;
      const beagle::discrete_dividend_schedule_t& dividends() const;
      const beagle::dividend_policy_ptr_t& dividendPolicy() const;
      const beagle::interp_builder_ptr_t& interpolation() const;
    private:
      beagle::dbl_t m_Spot;
      beagle::dbl_t m_Rate;
      beagle::dbl_t m_Volatility;
      int m_StepsPerAnnum;
      int m_StepsLogSpot;
      beagle::dbl_t m_NumStdev;
      beagle::discrete_dividend_schedule_t m_Dividends;
      beagle::dividend_policy_ptr_t m_Policy;
      beagle::interp_builder_ptr_t m_Interp;
    };

    struct Pricer
    {
      virtual ~Pricer( void );
    public:
      virtual beagle::dbl_t value( const beagle::product_ptr_t& product ) const = 0;
    public:
      static beagle::pricer_ptr_t formBlackScholesClosedFormEuropeanOptionPricer(
                                                                          beagle::dbl_t spot,
                                                                          beagle::dbl_t rate,
                                                                          beagle::dbl_t volatility,
                                                                          const discrete_dividend_schedule_t& dividends );
      static beagle::pricer_ptr_t formOneDimensionalBackwardPDEOptionPricer(
                                                                     const FiniteDifferenceDetails& fdDetails,
                                                                     const beagle::real_2d_function_ptr_t& diffusion );
      static beagle::pricer_ptr_t formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                     const FiniteDifferenceDetails& fdDetails,
                                                                     const beagle::real_2d_function_ptr_t& volatility);
      static beagle::pricer_ptr_t formOneDimForwardPDEEuroOptionPricer(beagle::dbl_t spot,
                                                                       const beagle::real_2d_function_ptr_t& drift,
                                                                       const beagle::real_2d_function_ptr_t& volatility,
                                                                       const beagle::real_function_ptr_t& rate,
                                                                       const beagle::discrete_dividend_schedule_t& dividends,
                                                                       const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
      static beagle::pricer_ptr_t formOneDimBackwardPDEOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                    const beagle::real_function_ptr_t& discounting,
                                                                    const beagle::real_2d_function_ptr_t& drift,
                                                                    const beagle::real_2d_function_ptr_t& volatility,
                                                                    const beagle::real_2d_function_ptr_t& rate,
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
        virtual void optionValueCollection( beagle::dbl_t start,
                                            beagle::dbl_t end,
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
