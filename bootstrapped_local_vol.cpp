#include "real_2d_function.hpp"
#include "payoff.hpp"
#include "pricer.hpp"
#include "interpolation_builder.hpp"
#include "real_function.hpp"

#include <iostream>

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct BootstrappedLocalVolatilityFunction : public RealTwoDimFunction
      {
        BootstrappedLocalVolatilityFunction( const beagle::dbl_vec_t& expiries,
                                             const beagle::dbl_vec_vec_t& strikesColl,
                                             const beagle::dbl_vec_vec_t& pricesColl,
                                             const beagle::pricer_ptr_t& forwardPricer,
                                             const beagle::payoff_ptr_t& payoff ) :
          m_Expiries( expiries ),
          m_StrikesColl( strikesColl ),
          m_PricesColl( pricesColl ),
          m_ForwardPricer( forwardPricer ),
          m_Payoff( payoff )
        { }
        virtual ~BootstrappedLocalVolatilityFunction( void )
        { }
      public:
        virtual double value( double argX,
                              double argY ) const override
        {
          if (!m_Func)
            doCalibration();

          return m_Func->value( argX, argY );
        }
      private:
        void doCalibration( void ) const
        {
          beagle::dbl_vec_t logStrikes;
          beagle::dbl_vec_t strikes;
          beagle::dbl_vec_t prices;
          beagle::real_function_ptr_coll_t localVols;

          auto pFD = dynamic_cast<beagle::valuation::mixins::FiniteDifference*>(m_ForwardPricer.get());
          double spot = pFD->spot();
          double rate = pFD->rate();
          int stepsPerAnnum = pFD->stepsPerAnnum();
          int numStrikes = pFD->numberOfStateVariableSteps();
          int numStdDev = pFD->numberOfStandardDeviations();
          const beagle::discrete_dividend_schedule_t& dividends = pFD->dividends();
          const beagle::dividend_policy_ptr_t& policy = pFD->dividendPolicy();
          const beagle::interp_builder_ptr_t& interp = pFD->interpolation();
          pFD->formStateVariableSteps( m_Expiries.back(), logStrikes, strikes );

          auto pOVCP = dynamic_cast<beagle::valuation::mixins::OptionValueCollectionProvider*>(m_ForwardPricer.get());
          if (!pOVCP)
            throw(std::string("Cast for OptionValueCollectionProvider failed!"));
          pOVCP->formInitialOptionValueCollection( m_Payoff, strikes, prices );

          double start(0.0);
          for (beagle::dbl_vec_t::size_type i=0; i<m_Expiries.size(); ++i)
          {
            double end = m_Expiries[i];
            beagle::dbl_vec_t initVols{.32, .29, .27, .25, .24, .235, .26};
            beagle::dbl_vec_t calibratedPrices(m_StrikesColl[i].size());
            beagle::real_function_ptr_t initGuess = interp->formFunction( m_StrikesColl[i], initVols );
            beagle::real_2d_function_ptr_t localVol 
                            = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( beagle::dbl_vec_t(1U, m_Expiries[i]),
                                                                                                      beagle::real_function_ptr_coll_t(1U, initGuess) );

            beagle::dbl_vec_t initPrices(strikes.size());
            do
            {
              std::copy( prices.cbegin(),
                         prices.cend(),
                         initPrices.begin() );
              beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                                        spot,
                                                                                        rate,
                                                                                        localVol,
                                                                                        stepsPerAnnum,
                                                                                        numStrikes,
                                                                                        numStdDev,
                                                                                        dividends,
                                                                                        policy,
                                                                                        interp );
              pOVCP = dynamic_cast<beagle::valuation::mixins::OptionValueCollectionProvider*>(forwardPricer.get());
              if (!pOVCP)
                throw(std::string("Cast for OptionValueCollectionProvider failed!"));

              pOVCP->optionValueCollection( start, end, m_Payoff, logStrikes, strikes, initPrices );
              beagle::real_function_ptr_t interpFunc = interp->formFunction( strikes, initPrices );

              for (beagle::dbl_vec_t::size_type j=0; j<m_StrikesColl[i].size(); ++j)
              {
                double calibratedPrice = interpFunc->value(m_StrikesColl[i][j]);
                calibratedPrices[j] = calibratedPrice;
                initVols[j] *= m_PricesColl[i][j] / calibratedPrice;
                // std::cout << m_PricesColl[i][j] << "\t" << calibratedPrice << std::endl;
              }
              // std::cout << std::endl;

              initGuess = interp->formFunction( m_StrikesColl[i], initVols );
              localVol = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( beagle::dbl_vec_t(1U, m_Expiries[i]),
                                                                                                 beagle::real_function_ptr_coll_t(1U, initGuess) );
            } while ( !hasConverged( m_PricesColl[i], calibratedPrices ) );

            for (beagle::dbl_vec_t::size_type j=0; j<calibratedPrices.size(); ++j)
            {
              std::cout << m_PricesColl[i][j] << "\t" << calibratedPrices[j] << std::endl;
            }
            std::cout << std::endl;

            start = end;
            prices = initPrices;
            localVols.push_back( initGuess );
          }

          m_Func = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( m_Expiries, localVols );
        }
      private:
        bool hasConverged( const beagle::dbl_vec_t& quotes,
                           const beagle::dbl_vec_t& prices ) const
        {
          double inner(0.0);
          for (beagle::dbl_vec_t::size_type i=0; i<quotes.size(); ++i)
          {
            inner += (quotes[i] - prices[i]) * (quotes[i] - prices[i]);
            std::cout << quotes[i] << "\t" << prices[i] << std::endl;
          }
          std::cout << std::endl;

          if (inner < 1e-3)
            return true;
          else
            return false;
        }
      private:
        beagle::dbl_vec_t m_Expiries;
        beagle::dbl_vec_vec_t m_StrikesColl;
        beagle::dbl_vec_vec_t m_PricesColl;
        beagle::pricer_ptr_t m_ForwardPricer;
        beagle::payoff_ptr_t m_Payoff;
        mutable beagle::real_2d_function_ptr_t m_Func;
      };

    }

    beagle::real_2d_function_ptr_t
    RealTwoDimFunction::createBootstrappedLocalVolatilityFunction( 
                                  const beagle::dbl_vec_t& expiries,
                                  const beagle::dbl_vec_vec_t& strikesColl,
                                  const beagle::dbl_vec_vec_t& pricesColl,
                                  const beagle::pricer_ptr_t& forwardPricer,
                                  const beagle::payoff_ptr_t& payoff )
    {
      if (expiries.size() != strikesColl.size())
        throw(std::string("Number of expiries must be identical to the number of strike collection"));
      if (expiries.size() != pricesColl.size())
        throw(std::string("Number of expiries must be identical to the number of price collection"));
      for (beagle::dbl_vec_vec_t::size_type i=0; i<strikesColl.size(); ++i)
      {
        if (strikesColl[i].size() != pricesColl[i].size())
          throw(std::string("The number of strikes must be identical to the number of prices"));
      }

      return std::make_shared<impl::BootstrappedLocalVolatilityFunction>( expiries,
                                                                          strikesColl,
                                                                          pricesColl,
                                                                          forwardPricer,
                                                                          payoff );
    }
  }
}