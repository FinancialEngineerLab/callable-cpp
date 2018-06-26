#include "euro_pricer.hpp"
#include "util.hpp"

namespace beagle
{
  namespace impl
  {
    BlackScholesClosedFormEuropeanOptionPricer::BlackScholesClosedFormEuropeanOptionPricer(
                                                    double spot,
                                                    double rate,
                                                    double volatility,
                                                    const discrete_dividend_schedule_t& dividends ) :
      m_Spot( spot ),
      m_Rate( rate ),
      m_Volatility( volatility ),
      m_Dividends( dividends )
    { }

    BlackScholesClosedFormEuropeanOptionPricer::~BlackScholesClosedFormEuropeanOptionPricer( void )
    { }

    void
    BlackScholesClosedFormEuropeanOptionPricer::calculateAdjustedSpotAndStrike( double expiry,
                                                                                double strike,
                                                                                double& adjustedSpot,
                                                                                double& adjustedStrike  ) const
    {
      adjustedSpot = m_Spot;
      adjustedStrike = strike;

      if (m_Dividends.empty())
        return;
      else
      {
        double expRateT = std::exp( m_Rate * expiry );
        double rootT = std::sqrt( expiry );
        double sigmaRootT = m_Volatility * rootT;
        double dOne = std::log( m_Spot * expRateT / strike ) / sigmaRootT + .5 * sigmaRootT;
        double dTwo = dOne - sigmaRootT;

        double nDOne = util::cumulativeStandardNormal(dOne);
        double nDTwo = util::cumulativeStandardNormal(dTwo);
        double denominator = nDOne - nDTwo;

        double phiDOne = util::standardNormal(dOne);
        double phiDTwo = util::standardNormal(dTwo);
        double phiDiff = phiDOne - phiDTwo;
        double crossDiff = phiDOne * nDTwo - phiDTwo * nDOne;

        double gamma = m_Spot * sigmaRootT * phiDOne * denominator * denominator * denominator;
        double a = - crossDiff * crossDiff;
        double b = phiDiff * crossDiff;
        double c = - phiDiff * phiDiff;
        double d = phiDOne * denominator * denominator;

        auto numDivs = m_Dividends.size();
        for (discrete_dividend_schedule_t::size_type i = 0U; i < numDivs; ++i)
        {
          double exDivTimeI = m_Dividends[i].first;
          double divAmountI = m_Dividends[i].second;

          if (exDivTimeI <= expiry)
          {
            /// First order term
            double discountingExDivTimeI = std::exp( - m_Rate * exDivTimeI );
            double dI = dOne - m_Volatility * exDivTimeI / rootT;
            double nDI = util::cumulativeStandardNormal(dI);

            double numeratorA =  nDI - nDTwo;
            double numeratorB = denominator - numeratorA;

            double aI = - discountingExDivTimeI * numeratorA / denominator;
            double bI = expRateT * discountingExDivTimeI * numeratorB / denominator;

            adjustedSpot += aI * divAmountI;
            adjustedStrike += bI * divAmountI;

            for (discrete_dividend_schedule_t::size_type j = i; j < numDivs; ++j)
            {
              double exDivTimeJ = m_Dividends[j].first;
              double divAmountJ = m_Dividends[j].second;

              if (exDivTimeJ <= expiry)
              {
                double discountingExDivTimeJ = std::exp( - m_Rate * exDivTimeJ );
                double dJ = dOne - m_Volatility * exDivTimeJ / rootT;
                double nDJ = util::cumulativeStandardNormal(dJ);

                double aIJ = a + b * (nDI + nDJ) + c * nDI * nDJ;

                double dIJ = dOne - m_Volatility * (exDivTimeI + exDivTimeJ) / rootT;
                double phiDIJ = util::standardNormal(dIJ);
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

    double
    BlackScholesClosedFormEuropeanOptionPricer::optionValue( const beagle::option_ptr_t& option ) const
    {
      auto pE = dynamic_cast<beagle::mixins::European*>( option.get() );
      if (!pE)
        throw(std::string("Cannot valuate an option with non-European exercise style in closed form!"));

      double expiry = option->expiry();
      double strike = option->strike();

      double adjustedSpot;
      double adjustedStrike;
      calculateAdjustedSpotAndStrike( expiry, strike, adjustedSpot, adjustedStrike );

      double discounting = std::exp( - m_Rate * expiry );
      double result = util::bsCall( adjustedStrike, adjustedSpot / discounting, expiry, m_Volatility ) * discounting;

      if (option->payoff()->isCall())
        return result;
      else
        return result - adjustedSpot + adjustedStrike * discounting;
    }
  }

  beagle::pricer_ptr_t
  Pricer::formBlackScholesClosedFormEuropeanOptionPricer( double spot,
                                                          double rate,
                                                          double volatility,
                                                          const beagle::discrete_dividend_schedule_t& dividends )
  {
    return std::make_shared<impl::BlackScholesClosedFormEuropeanOptionPricer>( spot, rate, volatility, dividends );
  }
}