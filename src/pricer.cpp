#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"

#include <iostream>

namespace beagle
{
  namespace valuation
  {
    Pricer::~Pricer( void )
    { }

    OneDimFiniteDifferenceSettings::OneDimFiniteDifferenceSettings( void ) :
      m_NumTimeSteps(1500),
      m_NumUnderlyingSteps(1501),
      m_NumStdev(7.0),
      m_Interp(beagle::math::InterpolationBuilder::linear()),
      m_Policy(beagle::valuation::DividendPolicy::liquidator())
    { }

    OneDimFiniteDifferenceSettings::OneDimFiniteDifferenceSettings( int numTimeSteps,
                                                                    int numStateVariableSteps,
                                                                    beagle::dbl_t numGaussianStandardDeviations,
                                                                    const beagle::interp_builder_ptr_t& interp,
                                                                    const beagle::dividend_policy_ptr_t& policy ) :
      m_NumTimeSteps(numTimeSteps),
      m_NumUnderlyingSteps(numStateVariableSteps),
      m_NumStdev(numGaussianStandardDeviations),
      m_Interp(interp),
      m_Policy(policy)
    { }

    OneDimFiniteDifferenceSettings::OneDimFiniteDifferenceSettings( int numTimeSteps,
                                                                    int numStateVariableSteps,
                                                                    beagle::dbl_t numGaussianStandardDeviations ) :
      m_NumTimeSteps(numTimeSteps),
      m_NumUnderlyingSteps(numStateVariableSteps),
      m_NumStdev(numGaussianStandardDeviations),
      m_Interp(beagle::math::InterpolationBuilder::linear()),
      m_Policy(beagle::valuation::DividendPolicy::liquidator())
    { }

    int OneDimFiniteDifferenceSettings::numberOfTimeSteps( void ) const
    {
      return m_NumTimeSteps;
    }

    int OneDimFiniteDifferenceSettings::numberOfStateVariableSteps( void ) const
    {
      return m_NumUnderlyingSteps;
    }

    beagle::dbl_t OneDimFiniteDifferenceSettings::numberOfGaussianStandardDeviations( void ) const
    {
      return m_NumStdev;
    }

    const beagle::dividend_policy_ptr_t& 
    OneDimFiniteDifferenceSettings::dividendPolicy( void ) const
    {
      return m_Policy;
    }

    const beagle::interp_builder_ptr_t& 
    OneDimFiniteDifferenceSettings::interpolationMethod( void ) const
    {
      return m_Interp;
    }

    FiniteDifferenceDetails::FiniteDifferenceDetails( void )
    { }

    FiniteDifferenceDetails::FiniteDifferenceDetails( beagle::dbl_t spot,
                                                      beagle::dbl_t rate,
                                                      beagle::dbl_t volatility,
                                                      int stepsPerAnnum,
                                                      int stepsLogSpot,
                                                      beagle::dbl_t numStdev,
                                                      const beagle::discrete_dividend_schedule_t& dividends,
                                                      const beagle::dividend_policy_ptr_t& policy,
                                                      const beagle::interp_builder_ptr_t& interp ) :
      m_Spot(spot),
      m_Rate(rate),
      m_Volatility(volatility),
      m_StepsPerAnnum(stepsPerAnnum),
      m_StepsLogSpot(stepsLogSpot),
      m_NumStdev(numStdev),
      m_Dividends(dividends),
      m_Policy(policy),
      m_Interp(interp)
    { }

    void
    FiniteDifferenceDetails::formTimeSteps( beagle::dbl_t start,
                                            beagle::dbl_t end,
                                            beagle::dbl_vec_t& times,
                                            beagle::int_vec_t& exDividendIndices ) const
    {
      exDividendIndices.clear();

      beagle::dbl_t expiry = end - start;
      int numSteps = static_cast<int>(std::floor(expiry * m_StepsPerAnnum));
      if (m_Dividends.empty())
      {
        times.resize(numSteps + 1);

        for (int i=0; i<numSteps+1; ++i)
          times[i] = start + i * expiry / numSteps;
      }
      else
      {
        times.reserve(numSteps + 1 + m_Dividends.size());

        auto it = m_Dividends.cbegin();
        if (it->first < start)
          ++it;

        auto itEnd = m_Dividends.cend();
        for (int i=0, j=0; i<numSteps+1; ++i, ++j)
        {
          beagle::dbl_t time = start + i * expiry / numSteps;
          if (it != itEnd)
          {
            if (it->first < time)
            {
              times.push_back(it->first);
              exDividendIndices.push_back(j);
              ++j;
              times.push_back(time);
              ++it;
            }
            else if (it->first == time)
            {
              times.push_back(it->first);
              exDividendIndices.push_back(j);
              ++it;
            }
          }

          times.push_back(time);
        }

        times.shrink_to_fit();
      }
    }

    void
    FiniteDifferenceDetails::formStateVariableSteps( beagle::dbl_t expiry,
                                                      beagle::dbl_vec_t& logStateVariables,
                                                      beagle::dbl_vec_t& stateVariables ) const
    {
      logStateVariables.resize(m_StepsLogSpot);
      stateVariables.resize(m_StepsLogSpot -2);

      beagle::dbl_t forward = m_Spot * std::exp(m_Rate * expiry);
      beagle::dbl_t logSpot = std::log( m_Spot ) + (m_Rate - .5 * m_Volatility * m_Volatility) * expiry;
      int mid = m_StepsLogSpot / 2;
      beagle::dbl_t logStrikestep = 2. * m_NumStdev * m_Volatility * std::sqrt(expiry) / m_StepsLogSpot;
      logStateVariables.resize(m_StepsLogSpot);
      for (int i=0; i<m_StepsLogSpot; ++i)
        logStateVariables[i] = logSpot + (i-mid)*logStrikestep;

      std::transform( logStateVariables.cbegin()+1,
                      logStateVariables.cend()-1,
                      stateVariables.begin(),
                      [](beagle::dbl_t arg) {return std::exp(arg);} );
    }

    beagle::dbl_t FiniteDifferenceDetails::spot() const
    {
      return m_Spot;
    }

    beagle::dbl_t FiniteDifferenceDetails::rate() const
    {
      return m_Rate;
    }

    beagle::dbl_t FiniteDifferenceDetails::volatility() const
    {
      return m_Volatility;
    }

    int FiniteDifferenceDetails::stepsPerAnnum() const
    {
      return m_StepsPerAnnum;
    }

    int FiniteDifferenceDetails::numberOfStateVariableSteps() const
    {
      return m_StepsLogSpot;
    }

    beagle::dbl_t FiniteDifferenceDetails::numberOfStandardDeviations() const
    {
      return m_NumStdev;
    }

    const beagle::discrete_dividend_schedule_t& FiniteDifferenceDetails::dividends() const
    {
      return m_Dividends;
    }

    const beagle::dividend_policy_ptr_t& FiniteDifferenceDetails::dividendPolicy() const
    {
      return m_Policy;
    }

    const beagle::interp_builder_ptr_t& FiniteDifferenceDetails::interpolation() const
    {
      return m_Interp;
    }

    namespace mixins
    {
      OptionValueCollectionProvider::~OptionValueCollectionProvider( void )
      { }

      FiniteDifference::~FiniteDifference( void )
      { }

      CloneWithNewLocalVolatilitySurface::~CloneWithNewLocalVolatilitySurface( void )
      { }

      CloneWithNewModelParameters::~CloneWithNewModelParameters( void )
      { }
    }
  }
}
