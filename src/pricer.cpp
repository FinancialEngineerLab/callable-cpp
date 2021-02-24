#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "real_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"
#include "util.hpp"
#include "option.hpp"

#include <iostream>

namespace beagle
{
  namespace valuation
  {
    OneDimFiniteDifferenceSettings::OneDimFiniteDifferenceSettings( void ) :
      m_NumTimeSteps(365),
      m_NumUnderlyingSteps(1001),
      m_NumStdev(4.5),
      m_Interp(beagle::math::InterpolationBuilder::linear())
    { }

    OneDimFiniteDifferenceSettings::OneDimFiniteDifferenceSettings( int numTimeSteps,
                                                                    int numStateVariableSteps,
                                                                    double numGaussianStandardDeviations,
                                                                    const beagle::interp_builder_ptr_t& interp ) :
      m_NumTimeSteps(numTimeSteps),
      m_NumUnderlyingSteps(numStateVariableSteps),
      m_NumStdev(numGaussianStandardDeviations),
      m_Interp(interp)
    { }

    OneDimFiniteDifferenceSettings::OneDimFiniteDifferenceSettings( int numTimeSteps,
                                                                    int numStateVariableSteps,
                                                                    double numGaussianStandardDeviations ) :
      m_NumTimeSteps(numTimeSteps),
      m_NumUnderlyingSteps(numStateVariableSteps),
      m_NumStdev(numGaussianStandardDeviations),
      m_Interp(beagle::math::InterpolationBuilder::linear())
    { }

    int OneDimFiniteDifferenceSettings::numberOfTimeSteps( void ) const
    {
      return m_NumTimeSteps;
    }

    int OneDimFiniteDifferenceSettings::numberOfStateVariableSteps( void ) const
    {
      return m_NumUnderlyingSteps;
    }

    double OneDimFiniteDifferenceSettings::numberOfGaussianStandardDeviations( void ) const
    {
      return m_NumStdev;
    }

    const beagle::interp_builder_ptr_t&
    OneDimFiniteDifferenceSettings::interpolationMethod( void ) const
    {
      return m_Interp;
    }

    ClosedFormEuropeanOptionPricer::ClosedFormEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                   const beagle::real_function_ptr_t& discounting) :
      m_Forward(forward),
      m_Discounting(discounting)
    { }

    double ClosedFormEuropeanOptionPricer::impliedBlackScholesVolatility(const beagle::product_ptr_t& product) const
    {
      auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
      if (!pE)
        throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

      auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
      if (!pO)
        throw(std::string("The incoming product is not an option!"));

      double expiry = pO->expiry();
      double strike = pO->strike();
      double price = value(product);
      return beagle::util::impliedBlackVolatility(price,
                                                  pO->strike(),
                                                  pO->expiry(),
                                                  pO->payoff(),
                                                  m_Forward->value(expiry),
                                                  m_Discounting->value(expiry));
    }

    const beagle::real_function_ptr_t& ClosedFormEuropeanOptionPricer::forwardCurve(void) const
    {
      return m_Forward;
    }
    const beagle::real_function_ptr_t& ClosedFormEuropeanOptionPricer::discountCurve(void) const
    {
      return m_Discounting;
    }

    namespace util
    {
      void checkSABRParameters(double alpha,
                                double beta,
                                double rho,
                                double nu,
                                bool isFreeBoundarySABR)
      {
        double betaMax = isFreeBoundarySABR ? .5 : 1.;
        const double& eps = beagle::util::epsilon();
        if (alpha < eps)
          throw("The alpha parameter of the SABR model must be positive!");
        if (beta < eps || beta - betaMax > eps)
          throw("The beta parameter of the SABR model must be between zero and " + std::to_string(betaMax) + "!");
        if (rho + 1. < eps || rho - 1. > eps)
          throw("The rho parameter of the SABR model must be between minus one and one!");
        if (nu < eps)
          throw("The nu parameter of the SABR model must be positive!");
      }

      void checkCEVParameters(double beta,
                              double sigma,
                              bool isFreeBoundaryCEV)
      {
        double betaMax = isFreeBoundaryCEV ? .5 : 1.;
        const double& eps = beagle::util::epsilon();
        if (beta < eps || beta - betaMax > eps)
          throw("The beta parameter of the CEV model must be between zero and " + std::to_string(betaMax) + "!");
        if (sigma < eps)
          throw("The sigma parameter of the CEV model must be positive!");
      }
    }
  }
}
