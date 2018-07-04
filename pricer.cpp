#include "pricer.hpp"

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
                                       int stepsPerAnnum,
                                       const beagle::discrete_dividend_schedule_t& dividends,
                                       beagle::dbl_vec_t& times,
                                       beagle::int_vec_t& exDividendIndices ) const
      {
        exDividendIndices.clear();

        double expiry = end - start;
        int numSteps = std::floor(expiry * stepsPerAnnum);
        if (dividends.empty())
        {
          times.resize(numSteps + 1);

          for (int i=0; i<numSteps+1; ++i)
            times[i] = start + i * expiry / numSteps;
        }
        else
        {
          times.reserve(numSteps + 1 + dividends.size());

          auto it = dividends.cbegin();
          if (it->first < start)
            ++it;

          auto itEnd = dividends.cend();
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
      FiniteDifference::formStateVariableSteps( double centralValue,
                                                double rate,
                                                double vol,
                                                double expiry,
                                                double numStdev,
                                                int numSteps,
                                                beagle::dbl_vec_t& logStateVariables,
                                                beagle::dbl_vec_t& stateVariables ) const
      {
        logStateVariables.resize(numSteps);
        stateVariables.resize(numSteps-2);

        double logSpot = std::log( centralValue ) + (rate - .5 * vol * vol) * expiry;
        int mid = numSteps / 2;
        double logStrikestep = 2. * numStdev * vol * std::sqrt(expiry) / numSteps;
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