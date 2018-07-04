#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "util.hpp"
#include "dividend_policy.hpp"
#include "real_function.hpp"
#include "interpolation_builder.hpp"

#include <fstream>
std::ofstream out("interpolation.txt");

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct OneDimensionalForwardPDEEuropeanOptionPricer : public Pricer,
                                                            public beagle::valuation::mixins::OptionValueCollectionProvider,
                                                            public beagle::valuation::mixins::FiniteDifference
      {
        using two_dbl_t = std::pair<double, double>;

        OneDimensionalForwardPDEEuropeanOptionPricer( double spot,
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
        ~OneDimensionalForwardPDEEuropeanOptionPricer( void )
        { }
      public:
        virtual void formInitialOptionValueCollection( const beagle::payoff_ptr_t& payoff,
                                                       const beagle::dbl_vec_t& strikes,
                                                       beagle::dbl_vec_t& prices ) const
        {
          prices.resize( strikes.size() );
          std::transform( strikes.cbegin(),
                          strikes.cend(),
                          prices.begin(),
                          [&payoff, this](double strike) {return payoff->intrinsicValue(m_Spot, strike);}  );
        }
        virtual void optionValueCollection( double start,
                                            double end,
                                            const beagle::payoff_ptr_t& payoff,
                                            const beagle::dbl_vec_t& logStrikes,
                                            const beagle::dbl_vec_t& strikes,
                                            beagle::dbl_vec_t& prices ) const override
        {
          beagle::dbl_vec_t times;
          beagle::int_vec_t exDividendIndices;
          formTimeSteps( start, end, m_StepsPerAnnum, m_Dividends, times, exDividendIndices );

          // for (auto price : prices)
          //   out << price << " ";
          // out << std::endl;

          two_dbl_t boundaryStrikes = std::make_pair( std::exp(logStrikes.front()),
                                                      std::exp(logStrikes.back()) );

          int timeSteps = times.size();
          double deltaX = logStrikes[1] - logStrikes[0];

          int strikeSize = strikes.size();
          beagle::dbl_vec_t diag(strikeSize);
          beagle::dbl_vec_t lower(strikeSize);
          beagle::dbl_vec_t upper(strikeSize);

          auto it = m_Dividends.cbegin();
          if (it->first < start)
            ++it;
          
          auto jt = exDividendIndices.cbegin();
          auto jtEnd = exDividendIndices.cend();
          for (int i=0; i<timeSteps-1; ++i)
          {
            double thisTime = times[i+1];
            double deltaT = thisTime - times[i];
            for (int j=0; j<strikeSize; ++j)
            {
              double vol = m_Volatility->value(thisTime, strikes[j]);
              double volOverDeltaX = vol / deltaX;
              double volOverDeltaXSquared = volOverDeltaX * volOverDeltaX;
              double mu = m_Rate + .5 * vol * vol;
              double muOverDeltaX = mu / deltaX;
              diag[j]  = 1. + deltaT * volOverDeltaXSquared;
              upper[j] =   deltaT * .5 * (muOverDeltaX - volOverDeltaXSquared);
              lower[j] = - deltaT * .5 * (muOverDeltaX + volOverDeltaXSquared);
            }

            two_dbl_t boundaryValues = boundaryCondition( payoff,
                                                          boundaryStrikes,
                                                          end - thisTime );
            prices[0]            -= deltaT * lower[0] * boundaryValues.first;
            prices[strikeSize-1] -= deltaT * upper[strikeSize-1] * boundaryValues.second;

            // for (auto price : prices)
            //   out << price << " ";
            // out << std::endl;

            beagle::util::tridiagonalSolve( prices, diag, upper, lower );

            /// Ex-dividend date
            if (jt != jtEnd && *jt == i)
            {
              beagle::real_function_ptr_t interpFunc = m_Interp->formFunction( strikes, prices );
              beagle::dbl_vec_t shiftedStrikes(strikes.size());

              double dividendAmount = it->second;
              std::transform( strikes.cbegin(),
                              strikes.cend(),
                              shiftedStrikes.begin(),
                              [this, dividendAmount](double strike) { 
                                return strike + m_Policy->dividendAmount(m_Spot, dividendAmount);
                              } );

              std::transform( shiftedStrikes.cbegin(),
                              shiftedStrikes.cend(),
                              prices.begin(),
                              [&interpFunc](double strike) { 
                                return interpFunc->value(strike);
                              } );

              ++jt;
              ++it;
            }
          }
        }
        virtual double optionValue( const beagle::option_ptr_t& option ) const override
        {
          if (auto pA = dynamic_cast<beagle::option::mixins::American*>(option.get()))
            throw(std::string("Cannot price an American option with forward PDE European option pricer!"));

          double expiry = option->expiry();
          double strike = option->strike();
          const beagle::payoff_ptr_t& payoff = option->payoff();

          double forward = m_Spot * std::exp(m_Rate * expiry);
          double atmVol = m_Volatility->value( expiry, forward );

          beagle::dbl_vec_t logStrikes;
          beagle::dbl_vec_t strikes;
          formStateVariableSteps( m_Spot, m_Rate, atmVol, expiry, m_NumStdev, m_StepsLogSpot, logStrikes, strikes );

          beagle::dbl_vec_t prices;
          formInitialOptionValueCollection( payoff, strikes, prices );
          optionValueCollection( 0., expiry, payoff, logStrikes, strikes, prices );

          beagle::real_function_ptr_t interpResult = m_Interp->formFunction( strikes, prices );
          return interpResult->value(strike);
        }
      private:
        two_dbl_t boundaryCondition( const beagle::payoff_ptr_t& payoff,
                                     const two_dbl_t& boundaryStrikes,
                                     double timeToExpiry ) const
        {
          double minStrike = boundaryStrikes.first;
          double maxStrike = boundaryStrikes.second;

          double discounting = std::exp(-m_Rate * timeToExpiry);
          return std::make_pair( payoff->intrinsicValue( m_Spot, minStrike * discounting ),
                                 payoff->intrinsicValue( m_Spot, maxStrike * discounting ) );
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
    Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer( double spot,
                                                              double rate,
                                                              const beagle::real_2d_function_ptr_t& volatility,
                                                              int stepsPerAnnum,
                                                              int stepsLogSpot,
                                                              double numStdev,
                                                              const beagle::discrete_dividend_schedule_t& dividends,
                                                              const beagle::dividend_policy_ptr_t& policy,
                                                              const beagle::interp_builder_ptr_t& interp )
    {
      return std::make_shared<impl::OneDimensionalForwardPDEEuropeanOptionPricer>( spot,
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