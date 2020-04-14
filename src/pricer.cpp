#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"

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

    namespace mixins
    {
      OptionValueCollectionProvider::~OptionValueCollectionProvider( void )
      { }

      CloneWithNewLocalVolatilitySurface::~CloneWithNewLocalVolatilitySurface( void )
      { }

      CloneWithNewModelParameters::~CloneWithNewModelParameters( void )
      { }
    }
  }
}
