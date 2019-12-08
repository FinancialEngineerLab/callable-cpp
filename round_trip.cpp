#include "round_trip.hpp"
#include "payoff.hpp"
#include "option.hpp"
#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"

namespace beagle
{
  namespace test
  {
    void generateEuropeanMarketQuotes( beagle::dbl_vec_t& expiries,
                                       beagle::dbl_vec_vec_t& strikesColl,
                                       beagle::dbl_vec_vec_t& pricesColl )
    {
      double spot = 100.;
      double rate = .00;
      beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();

      expiries.clear();
      strikesColl.clear();
      pricesColl.clear();

      expiries.push_back(.3);

      beagle::dbl_vec_t strikes{90., 92.5, 95., 97.5, 100., 102.5, 105., 107.5, 110.};
      beagle::dbl_vec_t vols{.44, .395, .355, .32, .29, .265, .28, .30, .33};
      strikesColl.push_back(strikes);

      auto it = strikes.cbegin();
      auto itEnd = strikes.cend();
      auto jt = vols.cbegin();
      beagle::dbl_vec_t prices;
      for ( ; it != itEnd; ++it, ++jt)
      {
        beagle::option_ptr_t amerOption = beagle::option::Option::createEuropeanOption( expiries[0],
                                                                                        *it,
                                                                                        payoff );
        beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimensionalBackwardPDEOptionPricer(
                                                               spot,
                                                               rate,
                                                               beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(*jt),
                                                               1501,
                                                               1901,
                                                               7.5,
                                                               beagle::discrete_dividend_schedule_t(),
                                                               beagle::valuation::DividendPolicy::liquidator(),
                                                               beagle::math::InterpolationBuilder::linear() );
        prices.push_back(odbpop->optionValue(amerOption));
      }

      pricesColl.push_back( prices );
    }
  }
}
