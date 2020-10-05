#include "pricer.hpp"
#include "option.hpp"
#include "bond.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "interpolation_builder.hpp"
#include "util.hpp"
#include "convertible_calibration.hpp"
#include "andersen_buffum.hpp"

#include <iostream>
#include <fstream>

namespace beagle
{
  namespace test
  {
    void generateAndersenBuffumFigureTwo( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;
      double sigma = .3;

      // Call option valuation with discounting, funding, and volatility
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);

      // Generate implied volatilities for a series of expiries
      beagle::dbl_vec_t expiries{.25, .5, 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 30.};

      // Panel A
      double p = 2.;
      beagle::dbl_vec_t cs{0., .03, .05, .1};

      std::ofstream outA(".\\figure\\fig_2_panel_A.txt");
      outA << "[";
      for (double expiry : expiries)
        outA << expiry << ", ";
      outA << "]\n";

      for (double c : cs)
      {
        outA << "[";
        for (double expiry : expiries)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = drift;

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
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outA << beagle::util::impliedBlackVolatility(odbpop->value(euroOption),
                                                      strike,
                                                      strike,
                                                      expiry,
                                                      discounting ) << ", ";
        }

        outA << "]\n";
      }

      // Panel B
      double c = .05;
      beagle::dbl_vec_t ps{0., .5, 2.};

      std::ofstream outB(".\\figure\\fig_2_panel_B.txt");
      outB << "[";
      for (double expiry : expiries)
        outB << expiry << ", ";
      outB << "]\n";

      for (double p : ps)
      {
        outB << "[";
        for (double expiry : expiries)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = drift;

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
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outB << beagle::util::impliedBlackVolatility(odbpop->value(euroOption),
                                                      strike,
                                                      strike,
                                                      expiry,
                                                      discounting ) << ", ";
        }

        outB << "]\n";
      }
    }

    void generateAndersenBuffumFigureThree( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;
      double sigma = .3;
      double p = 2.;

      // Call option valuation with discounting, funding, and volatility
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);

      // Generate implied volatility smile
      beagle::dbl_vec_t cs{.02, .05, .1};
      beagle::dbl_vec_t moneynesses{.6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.5};

      // Panel A
      double expiry = .5;

      std::ofstream outA(".\\figure\\fig_3_panel_A.txt");
      outA << "[";
      for (double moneyness : moneynesses)
        outA << moneyness << ", ";
      outA << "]\n";

      for (double c : cs)
      {
        outA << "[";
        for (double moneyness : moneynesses)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = drift;

          double strike = forward->value(expiry) * moneyness;
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
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outA << beagle::util::impliedBlackVolatility(odbpop->value(euroOption),
                                                      strike,
                                                      forward->value(expiry),
                                                      expiry,
                                                      discounting ) << ", ";
        }

        outA << "]\n";
      }

      // Panel B
      expiry = 5.;

      std::ofstream outB(".\\figure\\fig_3_panel_B.txt");
      outB << "[";
      for (double moneyness : moneynesses)
        outB << moneyness << ", ";
      outB << "]\n";

      for (double c : cs)
      {
        outB << "[";
        for (double moneyness : moneynesses)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = drift;

          double strike = forward->value(expiry) * moneyness;
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
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outB << beagle::util::impliedBlackVolatility(odbpop->value(euroOption),
                                                      strike,
                                                      forward->value(expiry),
                                                      expiry,
                                                      discounting ) << ", ";
        }

        outB << "]\n";
      }
    }

    void generateAndersenBuffumFigureFour( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;
      double sigma = .3;
      double expiry = 1;
      double c = .05;

      // Call option valuation with discounting, funding, and volatility
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);

      // Generate implied volatility smile
      beagle::dbl_vec_t ps{0., 2.};
      beagle::dbl_vec_t moneynesses{.6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.5};

      // Panel A

      std::ofstream outA(".\\figure\\fig_4.txt");
      outA << "[";
      for (double moneyness : moneynesses)
        outA << moneyness << ", ";
      outA << "]\n";

      for (double p : ps)
      {
        outA << "[";
        for (double moneyness : moneynesses)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = drift;

          double strike = forward->value(expiry) * moneyness;
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
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outA << beagle::util::impliedBlackVolatility(odbpop->value(euroOption),
                                                      strike,
                                                      forward->value(expiry),
                                                      expiry,
                                                      discounting ) << ", ";
        }

        outA << "]\n";
      }
    }

    void generateAndersenBuffumFigureFive( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;
      double sigma = .3;

      // Call option valuation with discounting, funding, and volatility
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::digitalCall();
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);

      // Generate implied volatilities for a series of expiries
      beagle::dbl_vec_t expiries{.25, .5, 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 30.};

      // Panel A
      double p = 2.;
      beagle::dbl_vec_t cs{.02, .05, .1};

      std::ofstream outA(".\\figure\\fig_5_panel_A.txt");
      outA << "[";
      for (double expiry : expiries)
        outA << expiry << ", ";
      outA << "]\n";

      for (double c : cs)
      {
        outA << "[";
        for (double expiry : expiries)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = drift;

          beagle::product_ptr_t riskyBond = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                  0.,
                                                                                                  payoff );

          beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                    forward,
                                                                    discounting,
                                                                    drift,
                                                                    volatility,
                                                                    rate,
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outA << -std::log(odbpop->value(riskyBond)) / expiry - r << ", ";
        }

        outA << "]\n";
      }

      // Panel B
      double c = .05;
      beagle::dbl_vec_t ps{0., 2., 3.};

      std::ofstream outB(".\\figure\\fig_5_panel_B.txt");
      outB << "[";
      for (double expiry : expiries)
        outB << expiry << ", ";
      outB << "]\n";

      for (double p : ps)
      {
        outB << "[";
        for (double expiry : expiries)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = drift;

          beagle::product_ptr_t riskyBond = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                  0.,
                                                                                                  payoff );

          beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                    forward,
                                                                    discounting,
                                                                    drift,
                                                                    volatility,
                                                                    rate,
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outB << -std::log(odbpop->value(riskyBond)) / expiry - r << ", ";
        }

        outB << "]\n";
      }
    }

    void generateAndersenBuffumFigureSix( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;
      double p = 2.;
      double c = .1;

      // Call option valuation with discounting, funding, and volatility
      beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::digitalCall();
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));
      beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                [=](double time, double price){ return c * std::pow(price / spot, -p); } );
      beagle::real_2d_function_ptr_t rate = drift;

      // Generate implied volatilities for a series of expiries
      beagle::dbl_vec_t expiries{.25, .5, 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 30.};

      // Panel A
      beagle::dbl_vec_t sigmas{.15, .3, .5};

      std::ofstream outA(".\\figure\\fig_6.txt");
      outA << "[";
      for (double expiry : expiries)
        outA << expiry << ", ";
      outA << "]\n";

      for (double sigma : sigmas)
      {
        beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);

        outA << "[";
        for (double expiry : expiries)
        {
          beagle::product_ptr_t riskyBond = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                  0.,
                                                                                                  payoff );

          beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                                    forward,
                                                                    discounting,
                                                                    drift,
                                                                    volatility,
                                                                    rate,
                                                                    beagle::valuation::OneDimFiniteDifferenceSettings(365, 750, 4.5) );

          outA << -std::log(odbpop->value(riskyBond)) / expiry - r << ", ";
        }

        outA << "]\n";
      }
    }

    void generateAndersenBuffumFigureSeven( void )
    {
      double spot = 50;
      double r = .04;
      double q = .02;

      double sigma = .4;
      double c = .05;

      beagle::valuation::OneDimFiniteDifferenceSettings settings(52, 150, 4.5);

      // Model parameters
      beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-r * arg);});
      beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                spot,
                                                beagle::math::RealFunction::createUnaryFunction(
                                                [=](double arg) { return std::exp(-(r - q) * arg);}));

      int numExpiries = 120U;
      beagle::dbl_vec_t expiries(numExpiries);
      for (int i=0; i<numExpiries; ++i)
        expiries[i] = (i+1) * 10. / numExpiries;

      beagle::andersen_buffum_param_t quotes(numExpiries, beagle::two_dbl_t{sigma, c});

      std::ofstream outA(".\\figure\\fig_7_panel_A.txt");
      outA << "[";
      for (double expiry : expiries)
        outA << expiry << ", ";
      outA << "]\n";

      std::ofstream outB(".\\figure\\fig_7_panel_B.txt");
      outB << "[";
      for (double expiry : expiries)
        outB << expiry << ", ";
      outB << "]\n";

      beagle::dbl_vec_t ps{0., .5, 1., 2.};
      for (double p : ps)
      {
        beagle::andersen_buffum_curve_pair_t curves =
                beagle::calibration::util::createCalibratedAndersenBuffumParameters(forward,
                                                                                    discounting,
                                                                                    settings,
                                                                                    p,
                                                                                    expiries,
                                                                                    quotes);

        outA << "[";
        outB << "[";
        for (int i=0; i<numExpiries; ++i)
        {
          double expiry = expiries[i];
          double fwd = forward->value(expiry);
          outA << curves.first->value(expiry, fwd) << ", ";
          outB << curves.second->value(expiry, spot) << ", ";
        }
        outA << "]\n";
        outB << "]\n";
      }
    }

    void generateAndersenBuffumTableOneCalibratedPrice( void )
    {
      // Convertible bond, Case A
      {
        std::cout << "\nConvertible bond, Case A\n";
        double spot = 50;
        double r = .04;
        double q = .02;
        double sigma = .4;

        double c = .03;
        double rec = 0.4;

        double expiry = 10.;
        double coupon = .03;
        int frequency = 2;

        beagle::valuation::OneDimFiniteDifferenceSettings settings(52, 150, 4.5);

        beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-r * arg);});
        beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                  spot,
                                                  beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-(r - q) * arg);}));

        // Create a fixed coupon bond: 5-year maturity, 1.5% coupon, semi-annual
        beagle::callable_schedule_t callSchedule;
        callSchedule.emplace_back(5., beagle::math::RealFunction::createConstantFunction(107.), 10.);

        beagle::puttable_schedule_t putSchedule;
        putSchedule.emplace_back(6., 100.);
        putSchedule.emplace_back(8., 100.);

        beagle::real_function_ptr_t conversionRatio = beagle::math::RealFunction::createConstantFunction(1.);
        beagle::product_ptr_t fcb = beagle::product::bond::Bond::createFixedCouponBond(expiry, coupon, frequency);
        beagle::product_ptr_t cb = beagle::product::bond::Bond::createConvertibleBond(fcb, conversionRatio, callSchedule, putSchedule);

        int numExpiries = static_cast<int>(expiry * 12);
        beagle::dbl_vec_t expiries(numExpiries);
        for (int i=0; i<numExpiries; ++i)
          expiries[i] = (i+1) * expiry / numExpiries;

        beagle::andersen_buffum_param_t quotes(numExpiries, beagle::two_dbl_t{sigma, c});

        beagle::dbl_vec_t ps{0., .5, 1., 2.};
        for (double p : ps)
        {
          beagle::andersen_buffum_curve_pair_t curves =
                  beagle::calibration::util::createCalibratedAndersenBuffumParameters(forward,
                                                                                      discounting,
                                                                                      settings,
                                                                                      p,
                                                                                      expiries,
                                                                                      quotes);

          beagle::real_2d_function_ptr_t drift = curves.second;
          beagle::real_2d_function_ptr_t volatility = curves.first;
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
                                                                  settings );

          std::cout << "The bond floor is:             " << odbpbp->value(fcb) << " for p = " << p << "\n";
          std::cout << "The convertible bond price is: " << odbpbp->value(cb) << " for p = " << p << "\n\n";
        }

        std::cout << "\n";
      }

      // Convertible bond, Case B
      {
        std::cout << "\nConvertible bond, Case B\n";
        double spot = 50;
        double r = .04;
        double q = .02;
        double sigma = .25;

        double c = .02;
        double rec = 0.4;

        double expiry = 5.;
        double coupon = .015;
        int frequency = 2;

        beagle::valuation::OneDimFiniteDifferenceSettings settings(52, 150, 4.5);

        beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-r * arg);});
        beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                  spot,
                                                  beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-(r - q) * arg);}));

        // Create a fixed coupon bond: 5-year maturity, 1.5% coupon, semi-annual
        beagle::callable_schedule_t callSchedule;
        callSchedule.emplace_back(3., beagle::math::RealFunction::createConstantFunction(107.), 5.);

        beagle::puttable_schedule_t putSchedule;
        putSchedule.emplace_back(4., 100.);

        beagle::real_function_ptr_t conversionRatio = beagle::math::RealFunction::createConstantFunction(1.);
        beagle::product_ptr_t fcb = beagle::product::bond::Bond::createFixedCouponBond(expiry, coupon, frequency);
        beagle::product_ptr_t cb = beagle::product::bond::Bond::createConvertibleBond(fcb, conversionRatio, callSchedule, putSchedule);

        int numExpiries = static_cast<int>(expiry * 12);
        beagle::dbl_vec_t expiries(numExpiries);
        for (int i=0; i<numExpiries; ++i)
          expiries[i] = (i+1) * expiry / numExpiries;

        beagle::andersen_buffum_param_t quotes(numExpiries, beagle::two_dbl_t{sigma, c});

        beagle::dbl_vec_t ps{0., .5, 1., 2.};
        for (double p : ps)
        {
          beagle::andersen_buffum_curve_pair_t curves =
                  beagle::calibration::util::createCalibratedAndersenBuffumParameters(forward,
                                                                                      discounting,
                                                                                      settings,
                                                                                      p,
                                                                                      expiries,
                                                                                      quotes);

          beagle::real_2d_function_ptr_t drift = curves.second;
          beagle::real_2d_function_ptr_t volatility = curves.first;
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
                                                                  settings );

          std::cout << "The bond floor is:             " << odbpbp->value(fcb) << " for p = " << p << "\n";
          std::cout << "The convertible bond price is: " << odbpbp->value(cb) << " for p = " << p << "\n\n";
        }

        std::cout << "\n";
      }
    }

    void generateAndersenBuffumTableOneNaivePrice(void)
    {
      std::cout << "\nStart of Test 7:\n\n";

      // Convertible bond, Case A
      {
        std::cout << "\nConvertible bond, Case A\n";
        double spot = 50;
        double r = .04;
        double q = .02;
        double sigma = .4;

        double c = .03;
        double rec = 0.4;

        double expiry = 10.;
        double coupon = .03;
        int frequency = 2;

        beagle::valuation::OneDimFiniteDifferenceSettings settings(104, 250, 4.5);

        beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-r*arg); });
        beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                  spot,
                                                  beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-(r-q)*arg); }));
        beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);

        // Create a fixed coupon bond: 5-year maturity, 1.5% coupon, semi-annual
        beagle::callable_schedule_t callSchedule;
        callSchedule.emplace_back(5., beagle::math::RealFunction::createConstantFunction(107.), 10.);

        beagle::puttable_schedule_t putSchedule;
        putSchedule.emplace_back(6., 100.);
        putSchedule.emplace_back(8., 100.);

        beagle::real_function_ptr_t conversionRatio = beagle::math::RealFunction::createConstantFunction(1.);
        beagle::product_ptr_t fcb = beagle::product::bond::Bond::createFixedCouponBond(expiry, coupon, frequency);
        beagle::product_ptr_t cb = beagle::product::bond::Bond::createConvertibleBond(fcb, conversionRatio, callSchedule, putSchedule);

        beagle::dbl_vec_t ps{0., .5, 1., 2.};
        for (double p : ps)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t recovery = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return -100. * rec * c * std::pow(price / spot, -p); } );
          beagle::pricer_ptr_t odbpbp  = beagle::valuation::Pricer::formOneDimBackwardPDEBondPricer(
                                                                  forward,
                                                                  discounting,
                                                                  drift,
                                                                  volatility,
                                                                  rate,
                                                                  recovery,
                                                                  settings );

          std::cout << "The bond floor is:             " << odbpbp->value(fcb) << " for p = " << p << "\n";
          std::cout << "The convertible bond price is: " << odbpbp->value(cb) << " for p = " << p << "\n\n";
        }

        std::cout << "\n";
      }

      // Convertible bond, Case B
      {
        std::cout << "\nConvertible bond, Case B\n";
        double spot = 50;
        double r = .04;
        double q = .02;
        double sigma = .25;

        double c = .02;
        double rec = 0.4;

        double expiry = 5.;
        double coupon = .015;
        int frequency = 2;

        beagle::valuation::OneDimFiniteDifferenceSettings settings(104, 250, 4.5);

        beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-r*arg); });
        beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                  spot,
                                                  beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-(r-q)*arg); }));
        beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);

        // Create a fixed coupon bond: 5-year maturity, 1.5% coupon, semi-annual
        beagle::callable_schedule_t callSchedule;
        callSchedule.emplace_back(3., beagle::math::RealFunction::createConstantFunction(105.), 5.);

        beagle::puttable_schedule_t putSchedule;
        putSchedule.emplace_back(4., 100.);

        beagle::real_function_ptr_t conversionRatio = beagle::math::RealFunction::createConstantFunction(1.);
        beagle::product_ptr_t fcb = beagle::product::bond::Bond::createFixedCouponBond(expiry, coupon, frequency);
        beagle::product_ptr_t cb = beagle::product::bond::Bond::createConvertibleBond(fcb, conversionRatio, callSchedule, putSchedule);

        beagle::dbl_vec_t ps{0., .5, 1., 2.};
        for (double p : ps)
        {
          beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t rate = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return c * std::pow(price / spot, -p); } );
          beagle::real_2d_function_ptr_t recovery = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                    [=](double time, double price){ return -100. * rec * c * std::pow(price / spot, -p); } );
          beagle::pricer_ptr_t odbpbp  = beagle::valuation::Pricer::formOneDimBackwardPDEBondPricer(
                                                                  forward,
                                                                  discounting,
                                                                  drift,
                                                                  volatility,
                                                                  rate,
                                                                  recovery,
                                                                  settings );

          std::cout << "The bond floor is:             " << odbpbp->value(fcb) << " for p = " << p << "\n";
          std::cout << "The convertible bond price is: " << odbpbp->value(cb) << " for p = " << p << "\n\n";
        }

        std::cout << "\n";
      }

      std::cout << "\nEnd of Test 7\n";
    }

    void test_volatility_smile_credit_spread_calibration(void)
    {
      std::cout << "\nStart of test calibrating Andersen-Buffum model to volatility smiles and credit spreads simutaneously:\n\n";

      // {
      //   // Forward is identical to spot
      //   double spot = 50;
      //   double r = .02;
      //   double q = .02;

      //   beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
      //                                             [=](double arg) { return std::exp(-r*arg); });
      //   beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
      //                                             spot,
      //                                             beagle::math::RealFunction::createUnaryFunction(
      //                                             [=](double arg) { return std::exp(-(r-q)*arg); }));
      //   beagle::valuation::OneDimFiniteDifferenceSettings settings(104, 250, 4.5);

      //   beagle::dbl_vec_t expiries{1./12, 2./12, 3./12, 6./12};
      //   beagle::dbl_vec_t strikes{42.5, 43.75, 45, 47.5, 50, 52.5, 55, 56.25, 57.5};
      //   beagle::dbl_vec_vec_t volatilities{{.347962, .333534, .319673, .294025, .272111, .255533, .245781, .243594, .243040},
      //                                      {.350014, .335503, .321562, .295765, .273724, .257051, .247243, .245043, .244487},
      //                                      {.352066, .337472, .323451, .297506, .275337, .258568, .248704, .246493, .245934},
      //                                      {.358224, .343380, .329118, .302727, .280177, .263120, .253089, .250842, .250276}};
      //   beagle::dbl_vec_t spreads{.03, .02, .01, .04};

      //   beagle::volatility_smile_credit_spread_coll_t quotes;
      //   for (beagle::dbl_vec_t::size_type i=0; i<expiries.size(); ++i)
      //   {
      //     double expiry = expiries[i];
      //     quotes.emplace_back(expiry,
      //                         std::make_pair(strikes, volatilities[i]),
      //                         spreads[i]);
      //   }

      //   beagle::dbl_vec_t ps{0., .5, 1., 2.};
      //   for (double p : ps)
      //   {
      //     beagle::andersen_buffum_curve_pair_t curves =
      //       beagle::calibration::util::createCalibratedAndersenBuffumParameters(forward,
      //                                                                           discounting,
      //                                                                           settings,
      //                                                                           p,
      //                                                                           quotes,
      //                                                                           beagle::math::InterpolationBuilder::piecewiseConstantRight());

      //     beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
      //                                                               forward,
      //                                                               discounting,
      //                                                               curves.second,
      //                                                               curves.first,
      //                                                               curves.second,
      //                                                               settings );

      //     std::cout << "\np = " << p << "\n\n";
      //     for (beagle::dbl_vec_t::size_type i=0; i<expiries.size(); ++i)
      //     {
      //       double expiry = expiries[i];
      //       double df = discounting->value(expiry);

      //       std::cout << "expiry = " << expiry;
      //       for (beagle::dbl_vec_t::size_type j=0; j<strikes.size(); ++j)
      //       {
      //         double strike = strikes[j];
      //         std::cout << "\n" << strike << "    " << curves.first->value(expiry, strike);
      //         beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
      //                                                                                                   strike,
      //                                                                                                   beagle::product::option::Payoff::call() );
      //         std::cout << "    " << odfpeop->value(euroOption)
      //                   << "    " << df * beagle::util::bsCall(strike, spot, expiry, volatilities[i][j]);
      //       }

      //       std::cout << "\n\n" << curves.second->value(expiry, spot);
      //       beagle::product_ptr_t riskyBond = beagle::product::option::Option::createEuropeanOption( expiry,
      //                                                                                                0.,
      //                                                                                                beagle::product::option::Payoff::digitalCall() );
      //       std::cout << "    " << odfpeop->value(riskyBond)
      //                 << "    " << df * std::exp(-spreads[i] * expiry) << "\n\n";
      //     }
      //   }
      // }

      // Conan's setup
      {
        double spot = 100;
        double r = .04;
        double q = -.01;

        beagle::real_function_ptr_t discounting = beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-r*arg); });
        beagle::real_function_ptr_t forward = beagle::math::RealFunction::createContinuousForwardAssetPriceFunction(
                                                  spot,
                                                  beagle::math::RealFunction::createUnaryFunction(
                                                  [=](double arg) { return std::exp(-(r-q)*arg); }));
        beagle::valuation::OneDimFiniteDifferenceSettings settings(104, 250, 4.5);

        beagle::dbl_vec_t expiries{1., 2.};
        beagle::dbl_vec_vec_t strikes{{105, 125}, {110, 128}};
        beagle::dbl_vec_vec_t volatilities{{.3753, .3253},
                                           {.4053, .3653}};
        beagle::dbl_vec_t spreads{.1, .1};

        beagle::volatility_smile_credit_spread_coll_t quotes;
        for (beagle::dbl_vec_t::size_type i=0; i<expiries.size(); ++i)
        {
          double expiry = expiries[i];
          quotes.emplace_back(expiry,
                              std::make_pair(strikes[i], volatilities[i]),
                              spreads[i]);
        }

        beagle::dbl_vec_t ps{1.5};
        for (double p : ps)
        {
          beagle::andersen_buffum_curve_pair_t curves =
            beagle::calibration::util::createCalibratedAndersenBuffumParameters(forward,
                                                                                discounting,
                                                                                settings,
                                                                                p,
                                                                                quotes,
                                                                                beagle::math::InterpolationBuilder::piecewiseConstantRight());

          beagle::pricer_ptr_t odfpeop  = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
                                                                    forward,
                                                                    discounting,
                                                                    curves.second,
                                                                    curves.first,
                                                                    curves.second,
                                                                    settings );

          std::cout << "\np = " << p << "\n\n";
          for (beagle::dbl_vec_t::size_type i=0; i<expiries.size(); ++i)
          {
            double expiry = expiries[i];
            double df = discounting->value(expiry);
            double fwd = forward->value(expiry);
            beagle::dbl_vec_t thisStrikes = strikes[i];

            std::cout << "expiry = " << expiry;
            for (beagle::dbl_vec_t::size_type j=0; j<thisStrikes.size(); ++j)
            {
              double strike = thisStrikes[j];
              std::cout << "\n" << strike << "    " << curves.first->value(expiry, strike);
              beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                        strike,
                                                                                                        beagle::product::option::Payoff::call() );
              std::cout << "    " << odfpeop->value(euroOption)
                        << "    " << df * beagle::util::bsCall(strike, fwd, expiry, volatilities[i][j]);
            }

            std::cout << "\n\n" << curves.second->value(expiry, spot);
            beagle::product_ptr_t riskyBond = beagle::product::option::Option::createEuropeanOption( expiry,
                                                                                                     0.,
                                                                                                     beagle::product::option::Payoff::digitalCall() );
            std::cout << "    " << odfpeop->value(riskyBond)
                      << "    " << df * std::exp(-spreads[i] * expiry) << "\n\n";
          }
        }
      }
      std::cout << "\nEnd of test\n";
    }
  }
}