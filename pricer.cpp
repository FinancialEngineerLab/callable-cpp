#include "pricer.hpp"
#include "real_2d_function.hpp"

namespace beagle
{
  namespace valuation
  {
    Pricer::Pricer( void )
    { }

    Pricer::~Pricer( void )
    { }

    namespace mixins
    {
      OptionValueCollectionProvider::~OptionValueCollectionProvider( void )
      { }

      FiniteDifference::~FiniteDifference( void )
      { }

      void
      FiniteDifference::formTimeSteps( double start,
                                       double end,
                                       beagle::dbl_vec_t& times,
                                       beagle::int_vec_t& exDividendIndices ) const
      {
        int steps = stepsPerAnnum();
        const beagle::discrete_dividend_schedule_t& divs = dividends();

        exDividendIndices.clear();

        double expiry = end - start;
        int numSteps = std::floor(expiry * steps);
        if (divs.empty())
        {
          times.resize(numSteps + 1);

          for (int i=0; i<numSteps+1; ++i)
            times[i] = start + i * expiry / numSteps;
        }
        else
        {
          times.reserve(numSteps + 1 + divs.size());

          auto it = divs.cbegin();
          if (it->first < start)
            ++it;

          auto itEnd = divs.cend();
          for (int i=0, j=0; i<numSteps+1; ++i, ++j)
          {
            double time = start + i * expiry / numSteps;
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
      FiniteDifference::formStateVariableSteps( double expiry,
                                                beagle::dbl_vec_t& logStateVariables,
                                                beagle::dbl_vec_t& stateVariables ) const
      {
        double centralValue = spot();
        double theRate = rate();
        const beagle::real_2d_function_ptr_t& vol = volatility();
        double numStdev = numberOfStandardDeviations();
        int numSteps = numberOfStateVariableSteps();

        logStateVariables.resize(numSteps);
        stateVariables.resize(numSteps-2);


        double forward = centralValue * std::exp(theRate * expiry);
        double atmVol = vol->value( expiry, forward );
        double logSpot = std::log( centralValue ) + (theRate - .5 * atmVol * atmVol) * expiry;
        int mid = numSteps / 2;
        double logStrikestep = 2. * numStdev * atmVol * std::sqrt(expiry) / numSteps;
        logStateVariables.resize(numSteps);
        for (int i=0; i<numSteps; ++i)
          logStateVariables[i] = logSpot + (i-mid)*logStrikestep;

        std::transform( logStateVariables.cbegin()+1,
                        logStateVariables.cend()-1,
                        stateVariables.begin(),
                        [](double arg) {return std::exp(arg);} );
      }
    }
  }
}