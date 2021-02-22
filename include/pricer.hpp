#pragma once

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
                                      double numGaussianStandardDeviations );
      OneDimFiniteDifferenceSettings( int numTimeSteps,
                                      int numStateVariableSteps,
                                      double numGaussianStandardDeviations,
                                      const beagle::interp_builder_ptr_t& interp );
    public:
      int numberOfTimeSteps( void ) const;
      int numberOfStateVariableSteps( void ) const;
      double numberOfGaussianStandardDeviations( void ) const;
      const beagle::interp_builder_ptr_t& interpolationMethod( void ) const;
    private:
      int m_NumTimeSteps;
      int m_NumUnderlyingSteps;
      double m_NumStdev;
      beagle::interp_builder_ptr_t m_Interp;
    };

    struct Pricer
    {
      virtual ~Pricer( void ) = default;
    public:
      virtual double value( const beagle::product_ptr_t& product ) const = 0;
    public:
      static beagle::pricer_ptr_t formBlackScholesClosedFormEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                                 const beagle::real_function_ptr_t& discounting,
                                                                                 double volatility );
      static beagle::pricer_ptr_t formOneDimForwardPDEArrowDebreuPricer(const beagle::real_function_ptr_t& forward,
                                                                        const beagle::real_function_ptr_t& discounting,
                                                                        const beagle::real_2d_function_ptr_t& drift,
                                                                        const beagle::real_2d_function_ptr_t& volatility,
                                                                        const beagle::real_2d_function_ptr_t& rate,
                                                                        const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
      static beagle::pricer_ptr_t formOneDimForwardPDEEuroOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                       const beagle::real_function_ptr_t& discounting,
                                                                       const beagle::real_2d_function_ptr_t& volatility,
                                                                       const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
      static beagle::pricer_ptr_t formOneDimBackwardPDEOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                    const beagle::real_function_ptr_t& discounting,
                                                                    const beagle::real_2d_function_ptr_t& drift,
                                                                    const beagle::real_2d_function_ptr_t& volatility,
                                                                    const beagle::real_2d_function_ptr_t& rate,
                                                                    const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
      static beagle::pricer_ptr_t formOneDimBackwardPDEBondPricer(const beagle::real_function_ptr_t& forward,
                                                                  const beagle::real_function_ptr_t& discounting,
                                                                  const beagle::real_2d_function_ptr_t& drift,
                                                                  const beagle::real_2d_function_ptr_t& volatility,
                                                                  const beagle::real_2d_function_ptr_t& rate,
                                                                  const beagle::real_2d_function_ptr_t& recovery,
                                                                  const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
      static beagle::pricer_ptr_t formClosedFormSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                         const beagle::real_function_ptr_t& discounting,
                                                                         const beagle::real_function_ptr_t& alpha,
                                                                         const beagle::real_function_ptr_t& beta,
                                                                         const beagle::real_function_ptr_t& rho,
                                                                         const beagle::real_function_ptr_t& nu );
      static beagle::pricer_ptr_t formClosedFormShiftedSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                                const beagle::real_function_ptr_t& discounting,
                                                                                const beagle::real_function_ptr_t& alpha,
                                                                                const beagle::real_function_ptr_t& beta,
                                                                                const beagle::real_function_ptr_t& rho,
                                                                                const beagle::real_function_ptr_t& nu,
                                                                                double shift );
      static beagle::pricer_ptr_t formClosedFormExactSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                              const beagle::real_function_ptr_t& discounting,
                                                                              const beagle::real_function_ptr_t& alpha,
                                                                              const beagle::real_function_ptr_t& beta,
                                                                              const beagle::real_function_ptr_t& rho,
                                                                              const beagle::real_function_ptr_t& nu,
                                                                              bool useApproximateKernel,
                                                                              const beagle::integration_method_ptr_t& quadMethod);
      static beagle::pricer_ptr_t formClosedFormFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                                     const beagle::real_function_ptr_t& discounting,
                                                                                     const beagle::real_function_ptr_t& alpha,
                                                                                     const beagle::real_function_ptr_t& beta,
                                                                                     const beagle::real_function_ptr_t& rho,
                                                                                     const beagle::real_function_ptr_t& nu,
                                                                                     bool useApproximateKernel,
                                                                                     const beagle::integration_method_ptr_t& quadMethod);
      static beagle::pricer_ptr_t formClsoedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                                                   const beagle::real_function_ptr_t& discounting,
                                                                                                   const beagle::real_function_ptr_t& alpha,
                                                                                                   const beagle::real_function_ptr_t& beta,
                                                                                                   const beagle::real_function_ptr_t& rho,
                                                                                                   const beagle::real_function_ptr_t& nu,
                                                                                                   bool useApproximateKernel,
                                                                                                   const beagle::integration_method_ptr_t& quadMethod );
      static beagle::pricer_ptr_t formClosedFormNormalFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                                           const beagle::real_function_ptr_t& discounting,
                                                                                           const beagle::real_function_ptr_t& alpha,
                                                                                           const beagle::real_function_ptr_t& rho,
                                                                                           const beagle::real_function_ptr_t& nu,
                                                                                           bool useApproximateKernel,
                                                                                           const beagle::integration_method_ptr_t& quadMethod );
      static beagle::pricer_ptr_t formClosedFormExactCEVEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                             const beagle::real_function_ptr_t& discounting,
                                                                             const beagle::real_function_ptr_t& beta,
                                                                             const beagle::real_function_ptr_t& sigma,
                                                                             const beagle::integration_method_ptr_t& quadMethod );
      static beagle::pricer_ptr_t formClosedFormFreeBoundaryCEVEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                                    const beagle::real_function_ptr_t& discounting,
                                                                                    const beagle::real_function_ptr_t& beta,
                                                                                    const beagle::real_function_ptr_t& sigma,
                                                                                    const beagle::integration_method_ptr_t& quadMethod);
    };

    namespace mixins
    {
      struct OneDimFokkerPlanck
      {
        virtual ~OneDimFokkerPlanck( void ) = default;
      public:
        virtual void formInitialCondition(double expiry,
                                          beagle::dbl_vec_t& stateVars,
                                          beagle::dbl_vec_t& density) const=0;
        virtual void evolve(double start,
                            double end,
                            const beagle::dbl_vec_t& stateVars,
                            beagle::dbl_vec_t& density) const=0;
      };

      struct Dupire
      {
        virtual ~Dupire( void ) = default;
      public:
        virtual void formInitialCondition(double expiry,
                                          const beagle::payoff_ptr_t& payoff,
                                          beagle::dbl_vec_t& stateVars,
                                          beagle::dbl_vec_t& prices) const=0;
        virtual void evolve(double start,
                            double end,
                            const beagle::payoff_ptr_t& payoff,
                            const beagle::dbl_vec_t& stateVars,
                            beagle::dbl_vec_t& prices) const=0;
      };

      struct ModelParameters
      {
        virtual ~ModelParameters(void) = default;
      public:
        virtual int numberOfParameters(void) const = 0;
        virtual double impliedBlackScholesVolatility(const beagle::product_ptr_t& product) const = 0;
      };
    }

    struct ClosedFormEuropeanOptionPricer : public Pricer,
                                            public mixins::ModelParameters
    {
    public:
      virtual int numberOfParameters(void) const = 0;
    };
  }
}
