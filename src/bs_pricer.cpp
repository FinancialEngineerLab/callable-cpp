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
          return beagle::util::bsValue(adjustedStrike, adjustedSpot * factor, expiry, m_Volatility, pO->payoff()) * discounting;
        }
      private:
        void calculateAdjustedSpotAndStrike( double expiry,
                                             double spot,
                                             double strike,
                                             double& adjustedSpot,
                                             double& adjustedStrike ) const
        {
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
            beagle::dbl_vec_t discreteDividends(diff);
            std::transform(dividends.cbegin(), it, discreteDividends.begin(),
                           [](const beagle::dividend_schedule_t::value_type& dividend)
                           { return std::get<2>(dividend); });

            bool zeroes = std::all_of(discreteDividends.cbegin(),
                                      discreteDividends.cend(),
                                      [](double dividend)
                                      { return std::fabs(dividend) < beagle::util::epsilon(); });
            if (zeroes)
              return;

            beagle::dbl_vec_t exDivTimes(diff);
            std::transform(dividends.cbegin(), it, exDivTimes.begin(),
                           [](const beagle::dividend_schedule_t::value_type& dividend)
                           { return std::get<0>(dividend); });

            auto pAF = dynamic_cast<beagle::math::mixins::AssetForward*>( m_Forward.get() );

            double discount = m_Discounting->value(expiry);
            double factor = pAF->multiplicativeForwardFactor(expiry);
            double rootT = std::sqrt( expiry );
            double sigmaRootT = m_Volatility * rootT;
            double dOne = std::log( spot * factor / strike ) / sigmaRootT + .5 * sigmaRootT;
            double dTwo = dOne - sigmaRootT;

            double nDOne = beagle::util::cumulativeStandardNormal(dOne);
            double nDTwo = beagle::util::cumulativeStandardNormal(dTwo);
            double denominator = nDOne - nDTwo;

            beagle::dbl_vec_t a(diff);
            beagle::dbl_vec_t b(diff);
            for (int i=0; i<diff; ++i)
            {
              double thisDiscount = m_Discounting->value(exDivTimes[i]);
              double thisFactor = pAF->multiplicativeForwardFactor(exDivTimes[i]);
              double dI = dOne - m_Volatility * exDivTimes[i] / rootT;
              double nDI = beagle::util::cumulativeStandardNormal(dI);

              a[i] = - (factor / thisFactor * nDI - thisDiscount / discount * nDTwo) / factor / denominator;
              b[i] =   (thisDiscount / discount * nDOne - factor / thisFactor * nDI) / denominator;

              adjustedSpot   += a[i] * discreteDividends[i];
              adjustedStrike += b[i] * discreteDividends[i];
            }

            double phiDOne = beagle::util::standardNormal(dOne);
            double phiDTwo = beagle::util::standardNormal(dTwo);
            double oneOverSpotVolRootT = 1. / spot / sigmaRootT;
            double prefactor = oneOverSpotVolRootT / denominator;
            for (int i=0; i<diff; ++i)
            {
              for (int j=i; j<diff; ++j)
              {
                double thisFactorI = pAF->multiplicativeForwardFactor(exDivTimes[i]);
                double thisFactorJ = pAF->multiplicativeForwardFactor(exDivTimes[j]);
                double dI = dOne - m_Volatility * (exDivTimes[i] + exDivTimes[j]) / rootT;
                double nDI = beagle::util::standardNormal(dI);

                double aIJ = 1. / thisFactorI / thisFactorJ * std::exp(m_Volatility*m_Volatility*exDivTimes[i]) * nDI
                           - (a[i] * phiDOne - b[i] / factor * phiDTwo) * (a[j] - spot * b[j] / strike);
                aIJ *= prefactor;
                double bIJ = aIJ * factor;

                if (i==j)
                {
                  aIJ /= 2.;
                  bIJ /= 2.;
                }

                adjustedSpot   += aIJ * discreteDividends[i] * discreteDividends[j];
                adjustedStrike += bIJ * discreteDividends[i] * discreteDividends[j];

              }
            }
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
