#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "util.hpp"
#include "dividend_policy.hpp"
#include "real_function.hpp"
#include "interpolation_builder.hpp"

#include <iostream>

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct OneDimensionalBackwardPDEOptionPricer : public Pricer
      {
        using two_dbl_t = std::pair<double, double>;

        OneDimensionalBackwardPDEOptionPricer( double spot,
                                               double rate,
                                               const beagle::real_2d_function_ptr_t& volatility,
                                               int stepsPerAnnum,
                                               int stepsLogSpot,
                                               double numStdev,
                                               const beagle::discrete_dividend_schedule_t& dividends,
                                               const beagle::dividend_policy_ptr_t& policy,
                                               const beagle::interp_builder_ptr_t& interp ) :
          m_Spot( spot ),
          m_Rate( rate ),
          m_Volatility( volatility ),
          m_StepsPerAnnum( stepsPerAnnum ),
          m_StepsLogSpot( stepsLogSpot ),
          m_NumStdev( numStdev ),
          m_Dividends( dividends ),
          m_Policy( policy ),
          m_Interp( interp )
        { }
        ~OneDimensionalBackwardPDEOptionPricer( void )
        { }
      public:
        virtual double optionValue( const beagle::option_ptr_t& option ) const override
        {
          double expiry = option->expiry();
          double strike = option->strike();
          const beagle::payoff_ptr_t& payoff = option->payoff();

          beagle::dbl_vec_t times;
          beagle::dbl_vec_t logSpots;
          beagle::int_vec_t exDividendIndices;
          formLatticeForBackwardValuation( expiry, times, logSpots, exDividendIndices );

          auto pA = dynamic_cast<beagle::option::mixins::American*>(option.get());
          bool isAmerican( pA != nullptr );

          int logSpotSize = logSpots.size();
          beagle::dbl_vec_t spots(logSpotSize);
          std::transform( logSpots.cbegin(),
                          logSpots.cend(),
                          spots.begin(),
                          [](double arg) {return std::exp(arg);} );

          // calculate terminal value for backward induction
          beagle::dbl_vec_t optionValues(logSpotSize);
          std::transform( spots.cbegin(),
                          spots.cend(),
                          optionValues.begin(),
                          [&payoff, strike](double spot) {return payoff->intrinsicValue(spot, strike);}  );

          int timeSteps = times.size();
          double deltaX = logSpots[1] - logSpots[0];

          beagle::dbl_vec_t diag(logSpotSize-2);
          beagle::dbl_vec_t lower(logSpotSize-2);
          beagle::dbl_vec_t upper(logSpotSize-2);
          beagle::dbl_vec_t rhs(optionValues.cbegin()+1, optionValues.cend()-1);

          auto it = m_Dividends.crbegin();
          auto jt = exDividendIndices.crbegin();
          auto jtEnd = exDividendIndices.crend();
          for (int i=timeSteps-1; i>0; --i)
          {
            double thisTime = times[i-1];
            double deltaT = times[i] - thisTime;
            for (int j=0; j<logSpotSize-2; ++j)
            {
              double vol = m_Volatility->value(thisTime, spots[j+1]);
              double volOverDeltaX = vol / deltaX;
              double volOverDeltaXSquared = volOverDeltaX * volOverDeltaX;
              double mu = m_Rate - .5 * vol * vol;
              double muOverDeltaX = mu / deltaX;
              diag[j]  = 1. + deltaT * (volOverDeltaXSquared + m_Rate);
              upper[j] = - deltaT * .5 * (muOverDeltaX + volOverDeltaXSquared);
              lower[j] =   deltaT * .5 * (muOverDeltaX - volOverDeltaXSquared);
            }

            two_dbl_t boundaryValues = boundaryCondition( payoff,
                                                          std::make_pair(spots.front(), spots.back()),
                                                          strike,
                                                          expiry - thisTime,
                                                          isAmerican );
            rhs[0]             -= deltaT * lower[0] * boundaryValues.first;
            rhs[logSpotSize-3] -= deltaT * upper[logSpotSize-3] * boundaryValues.second;

            beagle::util::tridiagonalSolve( rhs, diag, upper, lower );

            if (isAmerican)
            {
              for (int j=0; j<logSpotSize-2; ++j)
              {
                rhs[j] = std::max( payoff->intrinsicValue( spots[j+1], strike ), rhs[j] );
              }
            }

            /// Ex-dividend date
            if (jt != jtEnd && *jt == i)
            {
              beagle::dbl_vec_t shiftedSpots(spots.cbegin()+1, spots.cend()-1);
              beagle::real_function_ptr_t interpFunc = m_Interp->formFunction( shiftedSpots, rhs );

              double dividendAmount = it->second;
              std::transform( shiftedSpots.cbegin(),
                              shiftedSpots.cend(),
                              shiftedSpots.begin(),
                              [this, dividendAmount](double spot) { 
                                return m_Policy->exDividendStockPrice(spot, dividendAmount);
                              } );

              std::transform( shiftedSpots.cbegin(),
                              shiftedSpots.cend(),
                              rhs.begin(),
                              [&interpFunc](double spot) { 
                                return interpFunc->value(spot);
                              } );

              ++jt;
              ++it;
            }
          }

          beagle::dbl_vec_t spotsForInterp(spots.cbegin()+1, spots.cend()-1);
          beagle::real_function_ptr_t interpResult = m_Interp->formFunction( spotsForInterp, rhs );
          return interpResult->value(m_Spot);
        }
      private:
        void formLatticeForBackwardValuation( double expiry,
                                              beagle::dbl_vec_t& times,
                                              beagle::dbl_vec_t& logSpots,
                                              beagle::int_vec_t& exDividendIndices ) const
        {
          exDividendIndices.clear();

          int numSteps = std::floor(expiry * m_StepsPerAnnum);
          if (m_Dividends.empty())
          {
            times.resize(numSteps + 1);

            for (int i=0; i<numSteps+1; ++i)
              times[i] = i * expiry / numSteps;
          }
          else
          {
            times.reserve(numSteps + 1 + m_Dividends.size());

            auto it = m_Dividends.cbegin();
            auto itEnd = m_Dividends.cend();
            for (int i=0, j=0; i<numSteps+1; ++i, ++j)
            {
              double time = i * expiry / numSteps;
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
                else
                {
                  times.push_back(time);
                }
              }
            }

            times.shrink_to_fit();
          }

          double forward = m_Spot * std::exp(m_Rate * expiry);
          double atmVol = m_Volatility->value( expiry, forward );
          double logSpot = std::log( m_Spot ); // + (m_Rate - .5 * atmVol * atmVol) * expiry;
          int mid = m_StepsLogSpot / 2;
          double logSpotStep = 2. * m_NumStdev * atmVol * std::sqrt(expiry) / m_StepsLogSpot;
          logSpots.resize(m_StepsLogSpot);
          for (int i=0; i<m_StepsLogSpot; ++i)
            logSpots[i] = logSpot + (i-mid)*logSpotStep;
        }
        two_dbl_t boundaryCondition( const beagle::payoff_ptr_t& payoff,
                                     const two_dbl_t& boundarySpots,
                                     double strike,
                                     double timeToExpiry,
                                     bool isAmerican ) const
        {
          double minSpot = boundarySpots.first;
          double maxSpot = boundarySpots.second;

          double adjustedStrike = strike * std::exp(-m_Rate * timeToExpiry);
          if (isAmerican && payoff->isPut())
            adjustedStrike = strike;

          return std::make_pair( payoff->intrinsicValue( minSpot, adjustedStrike ),
                                 payoff->intrinsicValue( maxSpot, adjustedStrike ) );
        }
      private:
        double m_Spot;
        double m_Rate;
        beagle::real_2d_function_ptr_t m_Volatility;
        int m_StepsPerAnnum;
        int m_StepsLogSpot;
        double m_NumStdev;
        beagle::discrete_dividend_schedule_t m_Dividends;
        beagle::dividend_policy_ptr_t m_Policy;
        beagle::interp_builder_ptr_t m_Interp;
      };
    }

    beagle::pricer_ptr_t
    Pricer::formOneDimensionalBackwardPDEOptionPricer( double spot,
                                                       double rate,
                                                       const beagle::real_2d_function_ptr_t& volatility,
                                                       int stepsPerAnnum,
                                                       int stepsLogSpot,
                                                       double numStdev,
                                                       const beagle::discrete_dividend_schedule_t& dividends,
                                                       const beagle::dividend_policy_ptr_t& policy,
                                                       const beagle::interp_builder_ptr_t& interp )
    {
      return std::make_shared<impl::OneDimensionalBackwardPDEOptionPricer>( spot,
                                                                            rate,
                                                                            volatility,
                                                                            stepsPerAnnum,
                                                                            stepsLogSpot,
                                                                            numStdev,
                                                                            dividends,
                                                                            policy,
                                                                            interp );
    }
  }
}