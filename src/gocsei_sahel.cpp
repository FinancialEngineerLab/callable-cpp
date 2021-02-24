#include "gocsei_sahel.hpp"
#include "pricer.hpp"
#include "option.hpp"
#include "bond.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "dividend_policy.hpp"
#include "interpolation_builder.hpp"
#include "util.hpp"
#include "local_vol_calibration.hpp"

#include <iostream>
#include <fstream>


namespace beagle
{
  namespace test
  {
    void test_gocsei_sahel(void)
    {
      std::cout << "\nStart of Test 1:\n\n";

      // Set up dividend
      beagle::dividend_schedule_t dividends;
      dividends.emplace_back( 121. / 365., 0.0, 1.0 );
      //dividends.emplace_back( 1.5, 0.0, 3.0 );
      //dividends.emplace_back( 2.5, 0.01, .0 );
      //dividends.emplace_back( 3.5, 0.0, 3.0 );
      //dividends.emplace_back( 4.5, 0.03, .0 );
      //dividends.emplace_back( 5.5, 0.0, 3.0 );
      //dividends.emplace_back( 6.5, 0.0, 3.0 );
      //dividends.emplace_back( 7.5, 0.0, 3.0 );
      //dividends.emplace_back( 8.5, 0.0, 3.0 );
      //dividends.emplace_back( 9.5, 0.0, 3.0 );

      beagle::discrete_dividend_schedule_t discDivs(dividends.size());
      std::transform(dividends.cbegin(),
                    dividends.cend(),
                    discDivs.begin(),
                    [=](const beagle::dividend_schedule_t::value_type& item)
                    { return std::make_pair(std::get<0>(item), std::get<2>(item)); });

      // Model parameters
      double spot = 10.;
      double rate = .01;
      double carry = .0;
      double vol = .2;

      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-rate * arg);});
      beagle::real_function_ptr_t funding = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(rate - carry) * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createGeneralForwardAssetPriceFunction(
                                                spot,
                                                funding,
                                                dividends,
                                                beagle::valuation::DividendPolicy::liquidator());
      beagle::valuation::OneDimFiniteDifferenceSettings settings(750, 1501, 4.5);

      // Set up options
      double expiry = 364. / 365.;
      double strike = 10.;
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::put();

      beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                strike,
                                                                                                payoff );
      beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
                                                                                                strike,
                                                                                                payoff );

      try
      {
        beagle::pricer_ptr_t bscfeop = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(forward,
                                                                                                                 discounting,
                                                                                                                 vol );
        beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                  forward,
                                                                  discounting,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.),
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0),
                                                                  settings );
        beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                  forward,
                                                                  discounting,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vol),
                                                                  settings );

        std::cout << "European option price (CF)   is: " << bscfeop->value( euroOption ) << std::endl;
        std::cout << "European option price (FD-B) is: " << odbpop->value( euroOption ) << std::endl;
        std::cout << "American option price (FD-B) is: " << odbpop->value( amerOption ) << std::endl;
        std::cout << "European option price (FD-F) is: " << odfpeop->value( euroOption ) << std::endl;

        std::cout << "\nEnd of Test 1\n\n";
      }
      catch (const std::string& what)
      {
        std::cout << "A valuation error has occurred -- " << what << std::endl;
      }
    }

    void test_local_vol_calibration( void )
    {
      std::cout << "\nStart of Test 2:\n\n";

      {
        beagle::dividend_schedule_t dividends;
        dividends.emplace_back(0.6, 0.0, 3.0);
        dividends.emplace_back(1.4, 0.0, 3.0);

        beagle::dbl_vec_t strikes{90., 92.5, 95., 97.5, 100., 102.5, 105., 107.5, 110.};
        beagle::dbl_vec_t vols{.33, .32, .31, .30, .29, .28, .28, .29, .30};

        beagle::volatility_smile_coll_t volSmiles;
        volSmiles.emplace_back(0.25, std::make_pair(strikes, vols));
        volSmiles.emplace_back(0.50, std::make_pair(strikes, vols));
        volSmiles.emplace_back(0.75, std::make_pair(strikes, vols));
        volSmiles.emplace_back(1.00, std::make_pair(strikes, vols));
        volSmiles.emplace_back(1.25, std::make_pair(strikes, vols));
        volSmiles.emplace_back(1.50, std::make_pair(strikes, vols));
        volSmiles.emplace_back(1.75, std::make_pair(strikes, vols));
        volSmiles.emplace_back(2.00, std::make_pair(strikes, vols));

        double r = 0.04;
        double q = 0.02;
        double spot = 100;
        beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-r * arg);});
        beagle::real_function_ptr_t funding = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-(r-q) * arg);});
        beagle::real_function_ptr_t forward = beagle::math::RealFunction::createGeneralForwardAssetPriceFunction(
                                                  spot,
                                                  funding,
                                                  dividends,
                                                  beagle::valuation::DividendPolicy::liquidator());
        beagle::valuation::OneDimFiniteDifferenceSettings settings(104, 750, 10.5);

        beagle::real_2d_function_ptr_t localVol =
          beagle::calibration::util::createCalibratedLocalVolatilitySurface(forward,
                                                                            discounting,
                                                                            settings,
                                                                            volSmiles);
        //// Generate data for plot
        //beagle::dbl_vec_t expiries{.25, .5, .75, 1., 1.25, 1.5, 1.75, 2.};
        //beagle::dbl_vec_t expiriesPlot{0.};
        //for (double expiry : expiries)
        //{
        //  expiriesPlot.push_back(expiry - .00001);
        //  expiriesPlot.push_back(expiry);
        //  expiriesPlot.push_back(expiry + .00001);
        //}

        //beagle::dbl_vec_t strikesPlot(21U);
        //for (int i=0; i<21; ++i)
        //  strikesPlot[i] = 90 + i;

        //std::ofstream out(".\\figure\\AndreasenHuge\\local_vol_artificial3.txt");
        //out << "[";
        //for (double expiry : expiriesPlot)
        //  out << expiry << ", ";
        //out << "]\n";

        //out << "[";
        //for (double strike : strikesPlot)
        //  out << strike << ", ";
        //out << "]\n";

        //out << "[";
        //for (double expiry : expiriesPlot)
        //{
        //  out << "[";
        //  for (double strike : strikesPlot)
        //    out << localVol->value(expiry, strike) << ", ";
        //  out << "],\n";
        //}
        //out << "]";
      }

      {
        double spot = 2772.70;

        beagle::dbl_vec_t expiries{0.025, 0.101, 0.197, 0.274, 0.523, 0.772, 1.769, 2.267, 2.784, 3.781, 4.778, 5.774};

        beagle::dbl_vec_vec_t strikesColl;
        strikesColl.emplace_back(beagle::dbl_vec_t{0.8613, 0.8796, 0.8979, 0.9163, 0.9346, 0.9529, 0.9712, 0.9896, 1.0079, 1.0262, 1.0445, 1.0629, 1.0812, 1.0995, 1.1178});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.8796, 0.8979, 0.9163, 0.9346, 0.9529, 0.9712, 0.9896, 1.0079, 1.0262, 1.0445, 1.0629, 1.0812, 1.0995, 1.1178});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.8796, 0.8979, 0.9163, 0.9346, 0.9529, 0.9712, 0.9896, 1.0079, 1.0262, 1.0445, 1.0629, 1.0812, 1.0995, 1.1178});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.7697, 0.8063, 0.8430, 0.8796, 0.9163, 0.9529, 0.9896, 1.0262, 1.0629, 1.0995, 1.1362, 1.1728, 1.2095, 1.2461});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.7697, 0.8063, 0.8430, 0.8796, 0.9163, 0.9529, 0.9896, 1.0262, 1.0629, 1.0995, 1.1362, 1.1728, 1.2095, 1.2461});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.7697, 0.8063, 0.8430, 0.8796, 0.9163, 0.9529, 0.9896, 1.0262, 1.0629, 1.0995, 1.1362, 1.1728, 1.2095, 1.2461});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.7697, 0.8063, 0.8430, 0.8796, 0.9163, 0.9529, 0.9896, 1.0262, 1.0629, 1.0995, 1.1362, 1.1728, 1.2095, 1.2461});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.8063, 0.8796, 0.9529, 1.0262, 1.0995, 1.1728});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.5131, 0.5864, 0.6597, 0.7330, 0.8063, 0.8796, 0.9529, 1.0262, 1.0995, 1.1728, 1.2461, 1.3194, 1.3927, 1.4660});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.5131, 0.5864, 0.6597, 0.7330, 0.8063, 0.8796, 1.0262, 1.0995, 1.1728, 1.2461, 1.3194, 1.3927, 1.4660});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.6597, 0.7330, 0.8063, 0.8796, 0.9529, 1.0262, 1.0995, 1.1728, 1.2461, 1.3194, 1.3927, 1.4660});
        strikesColl.emplace_back(beagle::dbl_vec_t{0.8063, 0.8796, 0.9529, 1.0262, 1.0995, 1.1728, 1.2461, 1.3194, 1.3927});

        beagle::dbl_vec_vec_t volsColl;
        volsColl.emplace_back(beagle::dbl_vec_t{0.3365, 0.3216, 0.3043, 0.2880, 0.2724, 0.2586, 0.2466, 0.2358, 0.2247, 0.2159, 0.2091, 0.2056, 0.2045, 0.2025, 0.1933});
        volsColl.emplace_back(beagle::dbl_vec_t{0.2906, 0.2797, 0.2690, 0.2590, 0.2488, 0.2390, 0.2300, 0.2213, 0.2140, 0.2076, 0.2024, 0.1982, 0.1959, 0.1929});
        volsColl.emplace_back(beagle::dbl_vec_t{0.2764, 0.2672, 0.2578, 0.2489, 0.2405, 0.2329, 0.2253, 0.2184, 0.2123, 0.2069, 0.2025, 0.1984, 0.1944, 0.1920});
        volsColl.emplace_back(beagle::dbl_vec_t{0.3262, 0.3058, 0.2887, 0.2717, 0.2557, 0.2407, 0.2269, 0.2142, 0.2039, 0.1962, 0.1902, 0.1885, 0.1867, 0.1871});
        volsColl.emplace_back(beagle::dbl_vec_t{0.3079, 0.2936, 0.2798, 0.2663, 0.2531, 0.2404, 0.2284, 0.2173, 0.2074, 0.1988, 0.1914, 0.1854, 0.1811, 0.1785});
        volsColl.emplace_back(beagle::dbl_vec_t{0.3001, 0.2876, 0.2750, 0.2637, 0.2519, 0.2411, 0.2299, 0.2198, 0.2104, 0.2022, 0.1950, 0.1888, 0.1839, 0.1793});
        volsColl.emplace_back(beagle::dbl_vec_t{0.2843, 0.2753, 0.2666, 0.2575, 0.2497, 0.2418, 0.2347, 0.2283, 0.2213, 0.2151, 0.2091, 0.2039, 0.1990, 0.1945});
        volsColl.emplace_back(beagle::dbl_vec_t{0.2713, 0.2555, 0.2410, 0.2275, 0.2161, 0.2058});
        volsColl.emplace_back(beagle::dbl_vec_t{0.3366, 0.3178, 0.3019, 0.2863, 0.2711, 0.2580, 0.2448, 0.2322, 0.2219, 0.2122, 0.2054, 0.1988, 0.1930, 0.1849});
        volsColl.emplace_back(beagle::dbl_vec_t{0.3291, 0.3129, 0.2976, 0.2848, 0.2711, 0.2585, 0.2384, 0.2269, 0.2186, 0.2103, 0.2054, 0.2002, 0.1964});
        volsColl.emplace_back(beagle::dbl_vec_t{0.2975, 0.2848, 0.2722, 0.2611, 0.2501, 0.2392, 0.2305, 0.2223, 0.2164, 0.2105, 0.2054, 0.2012});
        volsColl.emplace_back(beagle::dbl_vec_t{0.2809, 0.2693, 0.2584, 0.2486, 0.2399, 0.2321, 0.2251, 0.2190, 0.2135});

        beagle::volatility_smile_coll_t volSmiles;
        for (beagle::dbl_vec_t::size_type i=0; i<expiries.size(); ++i)
        {
          beagle::dbl_vec_t strikes = strikesColl[i];
          for (beagle::dbl_vec_t::size_type j=0; j<strikes.size(); ++j)
            strikes[j] *= spot;

          volSmiles.emplace_back(expiries[i], std::make_pair(strikes, volsColl[i]));
        }

        double r = 0.0;
        double q = 0.0;
        beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-r * arg);});
        beagle::real_function_ptr_t funding = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-(r-q) * arg);});
        beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                  spot,
                                                  funding);
        beagle::valuation::OneDimFiniteDifferenceSettings settings(104, 750, 10.5);

        beagle::real_2d_function_ptr_t localVol =
          beagle::calibration::util::createCalibratedLocalVolatilitySurface(forward,
                                                                            discounting,
                                                                            settings,
                                                                            volSmiles);

        //// Generate data for plot
        //beagle::dbl_vec_t expiriesPlot{0.};
        //for (double expiry : expiries)
        //{
        //  expiriesPlot.push_back(expiry - .00001);
        //  expiriesPlot.push_back(expiry);
        //  expiriesPlot.push_back(expiry + .00001);
        //}
        //
        //beagle::dbl_vec_t strikesPlot(713U);
        //for (int i=0; i<713; ++i)
        //  strikesPlot[i] = 2388 + i;

        //std::ofstream out(".\\figure\\AndreasenHuge\\local_vol.txt");
        //out << "[";
        //for (double expiry : expiriesPlot)
        //  out << expiry << ", ";
        //out << "]\n";

        //out << "[";
        //for (double strike : strikesPlot)
        //  out << strike << ", ";
        //out << "]\n";

        //out << "[";
        //for (double expiry : expiriesPlot)
        //{
        //  out << "[";
        //  for (double strike : strikesPlot)
        //    out << localVol->value(expiry, strike) << ", ";
        //  out << "],\n";
        //}
        //out << "]";
      }

      std::cout << "\nEnd of Test 2\n\n";
    }

    void test_cev_finite_difference( void )
    {
      std::cout << "\nStart of Test 3:\n\n";

      beagle::discrete_dividend_schedule_t dividends;
      // dividends.emplace_back( 0.5, 6.0 );
      // dividends.emplace_back( 1.5, 6.5 );
      // dividends.emplace_back( 2.5, 7.0 );
      // dividends.emplace_back( 3.5, 7.5 );
      // dividends.emplace_back( 4.5, 8.0 );
      // dividends.emplace_back( 5.5, 8.0 );
      // dividends.emplace_back( 6.5, 8.0 );

      double expiry = 1.;
      double strike = 100.;
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::put();
      beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                strike,
                                                                                                payoff );
      beagle::product_ptr_t amerOption = beagle::product::option::Option::createAmericanOption( expiry,
                                                                                                strike,
                                                                                                payoff );

      double spot = 100.;
      double rate = .01;

      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return 1;});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot, discounting);

      // CEV modeling
      double alpha = .3;
      double beta = .5;
      beagle::real_function_ptr_t localVolFunction = beagle::math::RealFunction::createUnaryFunction(
                                                                  [alpha, beta](double x) { return alpha * std::pow(x, beta - 1.); } );
      beagle::real_2d_function_ptr_t localVolSurface =
                beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( beagle::dbl_vec_t(1U, 1.),
                                                                                        beagle::real_function_ptr_coll_t(1U, localVolFunction) );

      try
      {
        beagle::real_2d_function_ptr_t cev = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return alpha * std::pow(price, beta - 1.); } );
        beagle::pricer_ptr_t odbpop2  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                  forward,
                                                                  discounting,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  cev,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(500, 1001, 7.5) );
        beagle::pricer_ptr_t odfpeop2 = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
                                                                  forward,
                                                                  discounting,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  cev,
                                                                  beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(rate),
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(500, 1001, 7.5) );

        std::cout << "European option price (FD-B) is: " << odbpop2->value( euroOption ) << std::endl;
        std::cout << "American option price (FD-B) is: " << odbpop2->value( amerOption ) << std::endl;
        std::cout << "European option price (FD-F) is: " << odfpeop2->value( euroOption ) << std::endl;

        std::cout << "\nEnd of Test 3\n\n";
      }
      catch (const std::string& what)
      {
        std::cout << "A valuation error has occurred -- " << what << std::endl;
      }
    }

    void test_tridiagonal_solver(void)
    {
      std::cout << "\nStart of Test 4:\n\n";

      beagle::dbl_vec_t diag{ -2.6, -2.6, -2.6, -2.6 };
      beagle::dbl_vec_t upper{ 1., 1., 1., 1. };
      beagle::dbl_vec_t lower{ 1., 1., 1., 1. };
      beagle::dbl_vec_t rhs{ -240., 0., 0., -150. };

      beagle::util::tridiagonalSolve(rhs, diag, upper, lower);
      for (auto result : rhs)
        std::cout << result << "\n";

      std::cout << "\nEnd of Test 4:\n\n";
    }

    void test_implied_vol_credit_spread( void )
    {
      std::cout << "\nStart of Test 5:\n\n";

      double spot = 50;
      double r = .04;
      double q = .02;
      double sigma = .4;

      double c = .03;
      double p = 2.;

      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();

      // Model parameters
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                [=](double time, double price){ return c * std::pow(price / spot, -p); } );
      beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);
      beagle::real_2d_function_ptr_t rate = drift;

      beagle::dbl_vec_t expiries{5.}; //{.25, .5, 1., 2., 3., 4., 5.};
      for (double expiry : expiries)
      {
        double strike = forward->value(expiry);
        beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                  strike,
                                                                                                  payoff );
        beagle::product_ptr_t riskyBond = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                0.,
                                                                                                beagle::product::option::Payoff::digitalCall() );

        beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                  forward,
                                                                  discounting,
                                                                  drift,
                                                                  volatility,
                                                                  rate,
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(100, 500, 4.5) );
        double price = odbpop->value( euroOption );
        double bondPrice = odbpop->value(riskyBond);
        std::cout << expiry << "\t"
                  << price << "\t"
                  << beagle::util::impliedBlackVolatility(price,
                                                          strike,
                                                          expiry,
                                                          payoff,
                                                          forward->value(expiry),
                                                          discounting->value(expiry)) << "\t"
                  << bondPrice << "\t"
                  << -std::log(bondPrice) / expiry - r << "\n";

        beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
                                                                  forward,
                                                                  discounting,
                                                                  drift,
                                                                  volatility,
                                                                  rate,
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(100, 500, 4.5) );
        price = odfpeop->value( euroOption );
        bondPrice = odfpeop->value(riskyBond);
        std::cout << expiry << "\t"
                  << price << "\t"
                  << beagle::util::impliedBlackVolatility(price,
                                                          strike,
                                                          expiry,
                                                          payoff,
                                                          forward->value(expiry),
                                                          discounting->value(expiry)) << "\t"
                  << bondPrice << "\t"
                  << -std::log(bondPrice) / expiry - r << "\n\n";

      }

      std::cout << "\nEnd of Test 5\n\n";
    }

    void test_bond_pricer( void )
    {
      std::cout << "\nStart of Test 6:\n\n";

      // Model parameters 1
      double spot = 50;
      double r = .04;
      double q = .02;
      double sigma = .4;

      double c = .03;
      double p = 2.;
      double rec = 0.;

      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                [=](double time, double price){ return c * std::pow(price / spot, -p); } );
      beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);
      beagle::real_2d_function_ptr_t rate = drift;
      beagle::real_2d_function_ptr_t recovery = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                [=](double time, double price){ return -100. * rec * rate->value(time, price); } );
      beagle::pricer_ptr_t odbpbp  = beagle::valuation::Pricer::formOneDimBackwardPDEBondPricer(
                                                                  forward,
                                                                  discounting,
                                                                  drift,
                                                                  volatility,
                                                                  rate,
                                                                  recovery,
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(100, 500, 4.5) );

      // Create a fixed coupon bond: 5-year maturity, 1.5% coupon, semi-annual
      beagle::product_ptr_t fcb = beagle::product::bond::Bond::createFixedCouponBond(5., .015, 2);
      std::cout << "The bond price is: " << odbpbp->value(fcb) << "\n";

      // Model parameters 2
      sigma = .4;
      c = .03;
      p = 2.;

      drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                [=](double time, double price){ return c * std::pow(price / spot, -p); } );
      volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);
      rate = drift;
      recovery = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                [=](double time, double price){ return -100. * rec * rate->value(time, price); } );
      odbpbp  = beagle::valuation::Pricer::formOneDimBackwardPDEBondPricer(
                                                                  forward,
                                                                  discounting,
                                                                  drift,
                                                                  volatility,
                                                                  rate,
                                                                  recovery,
                                                                  beagle::valuation::OneDimFiniteDifferenceSettings(100, 500, 4.5) );

      // Create a fixed coupon bond: 10-year maturity, 3% coupon, semi-annual
      fcb = beagle::product::bond::Bond::createFixedCouponBond(10, .03, 2);
      std::cout << "The bond price is: " << odbpbp->value(fcb) << "\n";

      // Create a zero coupon bond: 5-year maturity
      beagle::product_ptr_t zcb = beagle::product::bond::Bond::createZeroCouponBond(5.);
      std::cout << "The bond price is: " << odbpbp->value(zcb) << "\n";

      std::cout << "\nEnd of Test 6\n\n";
    }

    void test_discontinuous_forward_curve( void )
    {
      std::cout << "\nStart of Test 8:\n\n";

      beagle::dividend_schedule_t dividends;
      dividends.emplace_back( 1.0, 0., 3.0 );
      dividends.emplace_back( 2.0, 0., 3.0 );
      dividends.emplace_back( 3.0, 0., 3.0 );
      dividends.emplace_back( 4.0, 0., 3.0 );
      dividends.emplace_back( 5.0, 0., 8.0 );
      dividends.emplace_back( 6.0, 0., 8.0 );

      double spot = 100.;
      double rate = .01;

      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-rate * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createGeneralForwardAssetPriceFunction(
                                                100,
                                                discounting,
                                                dividends,
                                                beagle::valuation::DividendPolicy::liquidator());

      beagle::dbl_vec_t expiries{.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
      for (double expiry : expiries)
      {
        std::cout << "Expiry = " << expiry << ", forward = " << forward->value(expiry) << "\n";
      }
      std::cout << "\n";

      std::cout << "\nEnd of Test 8\n\n";
    }

    void test_natural_cubic_spline( void )
    {
      std::cout << "\nStart of Test 9:\n\n";

      beagle::dbl_vec_t xValues{1., 2., 3., 4., 5.};
      beagle::dbl_vec_t yValues{0., 1., 0., 1., 0.};

      beagle::interp_builder_ptr_t spline = beagle::math::InterpolationBuilder::naturalCubicSpline();
      beagle::real_function_ptr_t func = spline->formFunction(xValues, yValues);

      beagle::dbl_vec_t xs(61U);
      std::cout << "[";
      for (int i=0; i<61; ++i)
      {
        xs[i] = i * .1;
        std::cout << xs[i] << ", ";
      }
      std::cout << "]\n";

      std::cout << "[";
      for (double x : xs)
      {
        std::cout << func->value(x) << ", ";
      }
      std::cout << "]\n";

      std::cout << "\nEnd of Test 9\n\n";
    }
  }
}