#include "pricer.hpp"
#include "util.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct BlackScholesClosedFormEuropeanOptionPricer : public Pricer
      {
        BlackScholesClosedFormEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                   const beagle::real_function_ptr_t& discounting,
                                                   double volatility) :
          m_Forward( forward ),
          m_Discounting( discounting ),
          m_Volatility( volatility )
        { }
        virtual ~BlackScholesClosedFormEuropeanOptionPricer( void ) = default;
      public:
        virtual double value( const beagle::product_ptr_t& product ) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();

          auto pAF = dynamic_cast<beagle::math::mixins::AssetForward*>( m_Forward.get() );
          if (!pAF)
            throw(std::string("Cannot find a forward curve!"));

          double spot = pAF->spot();
          double factor = pAF->multiplicativeForwardFactor(expiry);

          double adjustedSpot(spot);
          double adjustedStrike(strike);
          calculateAdjustedSpotAndStrike( expiry, spot, strike, adjustedSpot, adjustedStrike );

          double discounting = m_Discounting->value( expiry );
          double result = util::bsCall( adjustedStrike, adjustedSpot * factor, expiry, m_Volatility ) * discounting;

          if (pO->payoff()->isCall())
            return result;
          else
            return result - (adjustedSpot * factor - adjustedStrike) * discounting;
        }
      private:
        void calculateAdjustedSpotAndStrike( double expiry,
                                             double spot,
                                             double strike,
                                             double& adjustedSpot,
                                             double& adjustedStrike ) const
        {
          adjustedSpot = spot;
          adjustedStrike = strike;

          if (auto pCD = dynamic_cast<beagle::math::mixins::ContainsDividends*>(m_Forward.get()))
          {
            const beagle::dividend_schedule_t& dividends = pCD->dividendSchedule();
            auto it = std::upper_bound(dividends.cbegin(),
                                       dividends.cend(),
                                       expiry,
                                       [](double value,
                                          const beagle::dividend_schedule_t::value_type& dividend)
                                       { return value < std::get<0>(dividend); });

            auto diff = std::distance(dividends.cbegin(), it);
            beagle::dbl_vec_t exDivTimes(diff);
            std::transform(dividends.cbegin(), it, exDivTimes.begin(),
                           [](const beagle::dividend_schedule_t::value_type& dividend)
                           { return std::get<0>(dividend); });
            
            beagle::dbl_vec_t discreteDividends(diff);
            std::transform(dividends.cbegin(), it, discreteDividends.begin(),
                           [](const beagle::dividend_schedule_t::value_type& dividend)
                           { return std::get<2>(dividend); });

            auto pAF = dynamic_cast<beagle::math::mixins::AssetForward*>( m_Forward.get() );

            double discount = m_Discounting->value(expiry);
            double factor = pAF->multiplicativeForwardFactor(expiry);
            double rootT = std::sqrt( expiry );
            double sigmaRootT = m_Volatility * rootT;
            double dOne = std::log( spot * factor / strike ) / sigmaRootT + .5 * sigmaRootT;
            double dTwo = dOne - sigmaRootT;

            double nDOne = util::cumulativeStandardNormal(dOne);
            double nDTwo = util::cumulativeStandardNormal(dTwo);
            double denominator = nDOne - nDTwo;

            beagle::dbl_vec_t a(diff);
            beagle::dbl_vec_t b(diff);
            for (int i=0; i<diff; ++i)
            {
              double thisDiscount = m_Discounting->value(exDivTimes[i]);
              double thisFactor = pAF->multiplicativeForwardFactor(exDivTimes[i]);
              double dI = dOne - m_Volatility * exDivTimes[i] / rootT;
              double nDI = util::cumulativeStandardNormal(dI);

              a[i] = - (factor / thisFactor * nDI - thisDiscount / discount * nDTwo) / factor / denominator;
              b[i] =   (thisDiscount / discount * nDOne - factor / thisFactor * nDI) / denominator;

              adjustedSpot   = spot   + a[i] * discreteDividends[i];
              adjustedStrike = strike + b[i] * discreteDividends[i];
            }

          //  double phiDOne = util::standardNormal(dOne);
          //  double phiDTwo = util::standardNormal(dTwo);
          //  double phiDiff = phiDOne - phiDTwo;
          //  double crossDiff = phiDOne * nDTwo - phiDTwo * nDOne;

          //  double gamma = m_Spot * sigmaRootT * phiDOne * denominator * denominator * denominator;
          //  double a = - crossDiff * crossDiff;
          //  double b = phiDiff * crossDiff;
          //  double c = - phiDiff * phiDiff;
          //  double d = phiDOne * denominator * denominator;

          //  auto numDivs = m_Dividends.size();
          //  for (discrete_dividend_schedule_t::size_type i = 0U; i < numDivs; ++i)
          //  {
          //    double exDivTimeI = m_Dividends[i].first;
          //    double divAmountI = m_Dividends[i].second;

          //    if (exDivTimeI <= expiry)
          //    {
          //      /// First order term
          //      double discountingExDivTimeI = std::exp( - m_Rate * exDivTimeI );
          //      double dI = dOne - m_Volatility * exDivTimeI / rootT;
          //      double nDI = util::cumulativeStandardNormal(dI);

          //      double numeratorA =  nDI - nDTwo;
          //      double numeratorB = denominator - numeratorA;

          //      double aI = - discountingExDivTimeI * numeratorA / denominator;
          //      double bI = expRateT * discountingExDivTimeI * numeratorB / denominator;

          //      adjustedSpot += aI * divAmountI;
          //      adjustedStrike += bI * divAmountI;

          //      for (discrete_dividend_schedule_t::size_type j = i; j < numDivs; ++j)
          //      {
          //        double exDivTimeJ = m_Dividends[j].first;
          //        double divAmountJ = m_Dividends[j].second;

          //        if (exDivTimeJ <= expiry)
          //        {
          //          double discountingExDivTimeJ = std::exp( - m_Rate * exDivTimeJ );
          //          double dJ = dOne - m_Volatility * exDivTimeJ / rootT;
          //          double nDJ = util::cumulativeStandardNormal(dJ);

          //          double aIJ = a + b * (nDI + nDJ) + c * nDI * nDJ;

          //          double dIJ = dOne - m_Volatility * (exDivTimeI + exDivTimeJ) / rootT;
          //          double phiDIJ = util::standardNormal(dIJ);
          //          aIJ += d * std::exp(m_Volatility*m_Volatility*exDivTimeI) * phiDIJ;

          //          aIJ /= gamma;
          //          aIJ *= std::exp( -m_Rate * (exDivTimeI + exDivTimeJ));

          //          if (i == j)
          //            aIJ /= 2.;

          //          adjustedSpot += aIJ * divAmountI * divAmountJ;
          //          adjustedStrike += expRateT * aIJ * divAmountI * divAmountJ;
          //        }
          //      }
          //    }
          //  }
          }
        }
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        double m_Volatility;
      };
    }

    beagle::pricer_ptr_t
    Pricer::formBlackScholesClosedFormEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                           const beagle::real_function_ptr_t& discounting,
                                                           double volatility)
    {
      return std::make_shared<impl::BlackScholesClosedFormEuropeanOptionPricer>( forward, discounting, volatility );
    }
  }
}
