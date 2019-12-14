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

      expiries.push_back(.25);
      expiries.push_back(.5);
      //expiries.push_back(.75);
      //expiries.push_back(1.);

      beagle::dbl_vec_t strikes{90., 92.5, 95., 97.5, 100., 102.5, 105., 107.5, 110.};
      strikesColl.resize(expiries.size(), strikes);

      beagle::dbl_vec_t vols{.35, .32, .31, .30, .29, .28, .28, .29, .31};
      beagle::real_2d_function_ptr_t localVol
        = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(
                  expiries,
                  beagle::real_function_ptr_coll_t(
                              expiries.size(),
                              beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(strikes,
                                                                                                                vols)));

      beagle::payoff_ptr_t payoff = beagle::option::Payoff::call();

      beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                                 fdDetails,
                                                                 localVol );
      for (dbl_vec_t::size_type i = 0; i < expiries.size(); ++i)
      {
        beagle::dbl_vec_t prices;
        for (dbl_vec_t::size_type j = 0; j < strikesColl[i].size(); ++j)
        {
          beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(fdDetails.spot(),
                                                                                                                   fdDetails.rate(),
                                                                                                                   vols[j],
                                                                                                                   fdDetails.dividends());

          beagle::option_ptr_t euroOption = beagle::option::Option::createEuropeanOption(expiries[i],
                                                                                         strikesColl[i][j],
                                                                                         payoff);
          prices.push_back(odfpeop->optionValue(euroOption));
        }

        pricesColl.push_back(prices);
      }
    }
  }
}
