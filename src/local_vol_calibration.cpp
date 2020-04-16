#include "local_vol_calibration.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "calibration_adapter.hpp"
#include "calibration_functor.hpp"
#include "util.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "pricer.hpp"
#include "interpolation_builder.hpp"
#include "unsupported/Eigen/NonLinearOptimization"

#include <iostream>

namespace beagle
{
  namespace calibration
  {
    namespace util
    {
      namespace
      {
        struct LocalVolatilityCalibrationAdapter : public CalibrationAdapter
        {
          LocalVolatilityCalibrationAdapter(const beagle::real_function_ptr_t& forward,
                                            const beagle::real_function_ptr_t& discounting,
                                            const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                            const beagle::payoff_ptr_t& payoff,
                                            double start,
                                            double end,
                                            const beagle::dbl_vec_t& stateVars,
                                            const beagle::dbl_vec_t& prices,
                                            const beagle::volatility_smile_t& smile) :
            m_Forward(forward),
            m_Discounting(discounting),
            m_Settings(settings),
            m_Payoff(payoff),
            m_Start(start),
            m_End(end),
            m_StateVars(stateVars),
            m_Prices(prices)
          {
            m_Strikes = smile.first;
            beagle::dbl_vec_t vols = smile.second;

            m_Targets.clear();
            for (beagle::dbl_vec_t::size_type i=0; i<m_Strikes.size(); ++i)
            {
              double strike = m_Strikes[i];
              beagle::product_ptr_t euroOption = beagle::product::option::Option::createEuropeanOption( m_End,
                                                                                                        strike,
                                                                                                        payoff );
              beagle::pricer_ptr_t odbpop  = beagle::valuation::Pricer::formOneDimBackwardPDEOptionPricer(
                                                            forward,
                                                            discounting,
                                                            beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.),
                                                            beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(vols[i]),
                                                            beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0),
                                                            beagle::valuation::OneDimFiniteDifferenceSettings() );
              m_Targets.emplace_back(odbpop->value(euroOption));
            }

            std::cout << "\nFor this expiry, the targets are: \n";
            for (beagle::dbl_vec_t::size_type i=0; i<m_Targets.size(); ++i)
              std::cout << m_Targets[i] << "\t";
            std::cout << "\n";
          }
          virtual ~LocalVolatilityCalibrationAdapter( void ) = default;
        public:
          dbl_vec_t values( const dbl_vec_t& parameters ) const override
          {
            beagle::real_function_ptr_t localVol
                   = beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(m_Strikes, parameters);
            beagle::real_2d_function_ptr_t volatility
                   = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                   [=](double time, double spot){ return localVol->value(spot); });

            beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                     m_Forward,
                                                                     m_Discounting,
                                                                     volatility,
                                                                     m_Settings);
            auto pD = dynamic_cast<beagle::valuation::mixins::Dupire*>(forwardPricer.get());

            beagle::dbl_vec_t prices(m_Prices);
            pD->evolve(m_Start, m_End, m_Payoff, m_StateVars, prices);

            double fwd = m_Forward->value(m_End);
            double df = m_Discounting->value(m_End);
            beagle::dbl_vec_t strikes(m_StateVars);
            std::transform(m_StateVars.cbegin(),
                           m_StateVars.cend(),
                           strikes.begin(),
                           [=](double logMoneyness)
                           { return fwd * std::exp(logMoneyness); });
            beagle::real_function_ptr_t func = m_Settings.interpolationMethod()->formFunction(strikes, prices);
            
            beagle::dbl_vec_t results;
            for (beagle::dbl_vec_t::size_type j=0; j<m_Strikes.size(); ++j)
            {
              double strike = m_Strikes[j];
              results.emplace_back(df * fwd * func->value(strike) - m_Targets[j]);
            }

            return results;
          }
          dbl_mat_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const override
          {
            beagle::dbl_vec_t::size_type numRows = m_Strikes.size();
            dbl_mat_t result(numRows);
            for (beagle::dbl_vec_t::size_type i=0; i<numRows; ++i)
            {
              result[i].resize(parameters.size());
            }

            double bump = 1e-4;
            for (beagle::dbl_vec_t::size_type j=0; j<parameters.size(); ++j)
            {
              dbl_vec_t forwardParams(parameters);
              dbl_vec_t backwardParams(parameters);

              forwardParams[j] += bump;
              backwardParams[j] -= bump;

              dbl_vec_t forwardResult = values(forwardParams);
              dbl_vec_t backwardResult = values(backwardParams);

              for (beagle::dbl_vec_t::size_type i=0; i<numRows; ++i)
              {
                result[i][j] = ( forwardResult[i] - backwardResult[i] ) * .5 / bump;
              }
            }
            
            //// output parameters
            //std::cout << "\nFor this iteration, the calibrated parameters are: \n";
            //for (beagle::dbl_vec_t::size_type i=0; i<parameters.size(); ++i)
            //  std::cout << parameters[i] << "\t";
            //std::cout << "\n";

            // output residuals
            beagle::dbl_vec_t results = values(parameters);
            double max = 0.;
            std::cout << "\nFor this iteration, the relative errors are: \n";
            for (beagle::dbl_vec_t::size_type i=0; i<results.size(); ++i)
            {
              if (std::fabs(results[i] / m_Targets[i]) > max)
                max = results[i] / m_Targets[i];
              std::cout << results[i] / m_Targets[i] << "\t";
            }
            std::cout << "\n";
            std::cout << "The max relative error is: " << max << "\n";

            return result;
          }
        private:
          beagle::real_function_ptr_t m_Forward;
          beagle::real_function_ptr_t m_Discounting;
          beagle::valuation::OneDimFiniteDifferenceSettings m_Settings;
          beagle::payoff_ptr_t m_Payoff;
          double m_Start;
          double m_End;
          beagle::dbl_vec_t m_StateVars;
          beagle::dbl_vec_t m_Prices;
          beagle::dbl_vec_t m_Strikes;
          beagle::dbl_vec_t m_Targets;
        };
      }

      beagle::real_2d_function_ptr_t
      createCalibratedLocalVolatilitySurface(const beagle::real_function_ptr_t& forward,
                                             const beagle::real_function_ptr_t& discounting,
                                             const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                             const beagle::volatility_smile_coll_t& volSmiles)
      {
        // Set up the spatial finite difference grid and the initial condition
        int centralIndex = volSmiles.back().second.second.size() / 2;
        beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(
                                                            volSmiles.back().second.second[centralIndex]);

        beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                 forward,
                                                                 discounting,
                                                                 volatility,
                                                                 settings);
        auto pD = dynamic_cast<beagle::valuation::mixins::Dupire*>(forwardPricer.get());

        beagle::dbl_vec_t stateVars;
        beagle::dbl_vec_t prices;
        beagle::payoff_ptr_t call(beagle::product::option::Payoff::call());
        pD->formInitialCondition(volSmiles.back().first, call, stateVars, prices);

        // Now perform calibration by bootstrapping
        beagle::dbl_vec_t expiries;
        beagle::real_function_ptr_coll_t localVolFuncs;
        
        beagle::real_function_ptr_t localVol
                  = beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(volSmiles.front().second.first,
                                                                                                      volSmiles.front().second.second);

        double start = 0.;
        for (beagle::dbl_vec_t::size_type i=0; i<volSmiles.size(); ++i)
        {
          double end = volSmiles[i].first;
          beagle::dbl_vec_t strikes = volSmiles[i].second.first;

          beagle::calibration_bound_constraint_coll_t constraints(strikes.size(),
                                                                  beagle::calibration::CalibrationBoundConstraint::twoSidedBoundCalibrationConstraint(0.0001, 2.));
          beagle::int_vec_t elimIndices(0U);

          beagle::dbl_vec_t guesses(strikes);
          std::transform(strikes.cbegin(),
                         strikes.cend(),
                         guesses.begin(),
                         [=](double strike) { return localVol->value(strike); });

          guesses= beagle::calibration::util::getTransformedParameters( guesses, constraints );
          Eigen::VectorXd calibParams(guesses.size());
          for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
            calibParams(i) = guesses[i];

          beagle::calibration_adapter_ptr_t adapter = std::make_shared<LocalVolatilityCalibrationAdapter>(forward,
                                                                                                          discounting,
                                                                                                          settings,
                                                                                                          call,
                                                                                                          start,
                                                                                                          end,
                                                                                                          stateVars,
                                                                                                          prices,
                                                                                                          volSmiles[i].second);

          beagle::calibration::CalibrationFunctor functor( beagle::dbl_vec_t(strikes.size(), 0.),
                                                           adapter,
                                                           calibParams,
                                                           constraints,
                                                           elimIndices );
          Eigen::LevenbergMarquardt<beagle::calibration::CalibrationFunctor> lm(functor);
          lm.parameters.xtol = 1.0e-4;
          lm.parameters.maxfev = 50;
          Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(calibParams);

          for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
            guesses[i] = calibParams(i);

          guesses = beagle::calibration::util::getOriginalParameters( guesses, constraints );

          // Evolve density
          localVol = beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(strikes, guesses);
          beagle::real_2d_function_ptr_t volatility
                   = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                   [=](double time, double spot){ return localVol->value(spot); });
          forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                   forward,
                                                                   discounting,
                                                                   volatility,
                                                                   settings);
          auto pD = dynamic_cast<beagle::valuation::mixins::Dupire*>(forwardPricer.get());
          pD->evolve(start, end, call, stateVars, prices);

          start = end;
          expiries.emplace_back(end);
          localVolFuncs.emplace_back(localVol);

          std::cout << "\nCalibrated parameters are: \n";
          for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
            std::cout << guesses[i] << "\n";

          std::cout << "\n";
        }

        return beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(expiries, localVolFuncs);
      }
    }
  }
}