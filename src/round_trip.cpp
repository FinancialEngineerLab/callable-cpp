#include "round_trip.hpp"
#include "payoff.hpp"
#include "option.hpp"
#include "pricer.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"

namespace beagle
{
  namespace test
  {
    void generateEuropeanMarketQuotes( const beagle::valuation::FiniteDifferenceDetails& fdDetails,
                                       beagle::dbl_vec_t& expiries,
                                       beagle::dbl_vec_vec_t& strikesColl,
                                       beagle::dbl_vec_vec_t& pricesColl )
    {
      expiries.clear();
      strikesColl.clear();
      pricesColl.clear();

      expiries.push_back(.5);

      beagle::dbl_vec_t strikes{90., 92.5, 95., 97.5, 100., 102.5, 105., 107.5, 110.};
      strikesColl.push_back(strikes);

      beagle::dbl_vec_t vols{.44, .395, .355, .32, .29, .265, .28, .30, .33};
      beagle::real_2d_function_ptr_t localVol
        = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(
                  expiries,
                  beagle::real_function_ptr_coll_t(
                              1U,
                              beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(strikes,
                                                                                                                vols)));

      beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();

      beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                 fdDetails,
                                                                 localVol );

      auto it = strikes.cbegin();
      auto itEnd = strikes.cend();
      beagle::dbl_vec_t prices;
      for ( ; it != itEnd; ++it)
      {
        beagle::option_ptr_t amerOption = beagle::option::Option::createEuropeanOption( expiries[0],
                                                                                        *it,
                                                                                        payoff );
        prices.push_back(odfpeop->optionValue(amerOption));
      }

      pricesColl.push_back( prices );
    }
  }
}
