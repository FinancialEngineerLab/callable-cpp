#include "pricer.hpp"
#include "util.hpp"
#include "option.hpp"
#include "payoff.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct BlackScholesClosedFormEuropeanOptionPricer : public Pricer,
                                                          public beagle::valuation::mixins::CloneWithNewModelParameters
      {
        BlackScholesClosedFormEuropeanOptionPricer( beagle::dbl_t spot,
                                                    beagle::dbl_t rate,
                                                    beagle::dbl_t volatility,
                                                    const beagle::discrete_dividend_schedule_t& dividends ) :
          m_Spot( spot ),
          m_Rate( rate ),
          m_Volatility( volatility ),
          m_Dividends( dividends )
        { }
        ~BlackScholesClosedFormEuropeanOptionPricer( void )
        { }
      public:
        virtual beagle::dbl_t value( const beagle::product_ptr_t& product ) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          beagle::dbl_t expiry = pO->expiry();
          beagle::dbl_t strike = pO->strike();

          beagle::dbl_t adjustedSpot;
          beagle::dbl_t adjustedStrike;
          calculateAdjustedSpotAndStrike( expiry, strike, adjustedSpot, adjustedStrike );

          beagle::dbl_t discounting = std::exp( - m_Rate * expiry );
          beagle::dbl_t result = util::bsCall( adjustedStrike, adjustedSpot / discounting, expiry, m_Volatility ) * discounting;

          if (pO->payoff()->isCall())
            return result;
          else
            return result - adjustedSpot + adjustedStrike * discounting;
        }
        virtual beagle::pricer_ptr_t createPricerWithNewModelParameters( const beagle::dbl_vec_t& parameters ) const override
        {
          if (parameters.size() != 1U)
            throw("Cannot create a new Black-Scholes pricer with more than one volatility");

          return beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(m_Spot, m_Rate, parameters.front(), m_Dividends);
        }
      private:
        void calculateAdjustedSpotAndStrike( beagle::dbl_t expiry,
                                             beagle::dbl_t strike,
                                             beagle::dbl_t& adjustedSpot,
                                             beagle::dbl_t& adjustedStrike ) const
        {
          adjustedSpot = m_Spot;
          adjustedStrike = strike;

          if (m_Dividends.empty())
            return;
          else
          {
            beagle::dbl_t expRateT = std::exp( m_Rate * expiry );
            beagle::dbl_t rootT = std::sqrt( expiry );
            beagle::dbl_t sigmaRootT = m_Volatility * rootT;
            beagle::dbl_t dOne = std::log( m_Spot * expRateT / strike ) / sigmaRootT + .5 * sigmaRootT;
            beagle::dbl_t dTwo = dOne - sigmaRootT;

            beagle::dbl_t nDOne = util::cumulativeStandardNormal(dOne);
            beagle::dbl_t nDTwo = util::cumulativeStandardNormal(dTwo);
            beagle::dbl_t denominator = nDOne - nDTwo;

            beagle::dbl_t phiDOne = util::standardNormal(dOne);
            beagle::dbl_t phiDTwo = util::standardNormal(dTwo);
            beagle::dbl_t phiDiff = phiDOne - phiDTwo;
            beagle::dbl_t crossDiff = phiDOne * nDTwo - phiDTwo * nDOne;

            beagle::dbl_t gamma = m_Spot * sigmaRootT * phiDOne * denominator * denominator * denominator;
            beagle::dbl_t a = - crossDiff * crossDiff;
            beagle::dbl_t b = phiDiff * crossDiff;
            beagle::dbl_t c = - phiDiff * phiDiff;
            beagle::dbl_t d = phiDOne * denominator * denominator;

            auto numDivs = m_Dividends.size();
            for (discrete_dividend_schedule_t::size_type i = 0U; i < numDivs; ++i)
            {
              beagle::dbl_t exDivTimeI = m_Dividends[i].first;
              beagle::dbl_t divAmountI = m_Dividends[i].second;

              if (exDivTimeI <= expiry)
              {
                /// First order term
                beagle::dbl_t discountingExDivTimeI = std::exp( - m_Rate * exDivTimeI );
                beagle::dbl_t dI = dOne - m_Volatility * exDivTimeI / rootT;
                beagle::dbl_t nDI = util::cumulativeStandardNormal(dI);

                beagle::dbl_t numeratorA =  nDI - nDTwo;
                beagle::dbl_t numeratorB = denominator - numeratorA;

                beagle::dbl_t aI = - discountingExDivTimeI * numeratorA / denominator;
                beagle::dbl_t bI = expRateT * discountingExDivTimeI * numeratorB / denominator;

                adjustedSpot += aI * divAmountI;
                adjustedStrike += bI * divAmountI;

                for (discrete_dividend_schedule_t::size_type j = i; j < numDivs; ++j)
                {
                  beagle::dbl_t exDivTimeJ = m_Dividends[j].first;
                  beagle::dbl_t divAmountJ = m_Dividends[j].second;

                  if (exDivTimeJ <= expiry)
                  {
                    beagle::dbl_t discountingExDivTimeJ = std::exp( - m_Rate * exDivTimeJ );
                    beagle::dbl_t dJ = dOne - m_Volatility * exDivTimeJ / rootT;
                    beagle::dbl_t nDJ = util::cumulativeStandardNormal(dJ);

                    beagle::dbl_t aIJ = a + b * (nDI + nDJ) + c * nDI * nDJ;

                    beagle::dbl_t dIJ = dOne - m_Volatility * (exDivTimeI + exDivTimeJ) / rootT;
                    beagle::dbl_t phiDIJ = util::standardNormal(dIJ);
                    aIJ += d * std::exp(m_Volatility*m_Volatility*exDivTimeI) * phiDIJ;

                    aIJ /= gamma;
                    aIJ *= std::exp( -m_Rate * (exDivTimeI + exDivTimeJ));

                    if (i == j)
                      aIJ /= 2.;

                    adjustedSpot += aIJ * divAmountI * divAmountJ;
                    adjustedStrike += expRateT * aIJ * divAmountI * divAmountJ;
                  }
                }
              }
            }
          }
        }
      private:
        beagle::dbl_t m_Spot;
        beagle::dbl_t m_Rate;
        beagle::dbl_t m_Volatility;
        beagle::discrete_dividend_schedule_t m_Dividends;
      };
    }

    beagle::pricer_ptr_t
    Pricer::formBlackScholesClosedFormEuropeanOptionPricer( beagle::dbl_t spot,
                                                            beagle::dbl_t rate,
                                                            beagle::dbl_t volatility,
                                                            const beagle::discrete_dividend_schedule_t& dividends )
    {
      return std::make_shared<impl::BlackScholesClosedFormEuropeanOptionPricer>( spot, rate, volatility, dividends );
    }
  }
}
