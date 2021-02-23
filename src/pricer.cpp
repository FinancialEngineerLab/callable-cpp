#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"
#include "util.hpp"

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
