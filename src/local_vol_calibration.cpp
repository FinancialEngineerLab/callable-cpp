#include "local_vol_calibration.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "calibration_adapter.hpp"
#include "calibration_functor.hpp"
#include "util.hpp"
#include "payoff.hpp"
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
                                           double start,
                                           double end,
                                           const beagle::dbl_vec_t& stateVars,
                                           const beagle::dbl_vec_t& density,
                                           const beagle::volatility_smile_t& smile) :
            m_Forward(forward),
            m_Discounting(discounting),
            m_Settings(settings),
            m_Start(start),
            m_End(end),
            m_StateVars(stateVars),
            m_Density(density)
          {
            double df = discounting->value(end);
            double fwd = forward->value(end);

            m_Strikes = smile.first;
            beagle::dbl_vec_t vols = smile.second;

            m_Targets.clear();
            for (beagle::dbl_vec_t::size_type i=0; i<m_Strikes.size(); ++i)
            {
              double strike = m_Strikes[i];
              double call = beagle::util::bsCall(m_Strikes[i], fwd, end, vols[i]);
              if (strike > fwd)
                m_Targets.emplace_back(df * call);
              else
                m_Targets.emplace_back(df * (call - fwd + strike));
            }
          }
          virtual ~LocalVolatilityCalibrationAdapter( void )
          { }
        public:
          dbl_vec_t values( const dbl_vec_t& parameters ) const override
          {
            //std::cout << "\nFor this iteration, the calibrated parameters are: \n";
            //for (beagle::dbl_vec_t::size_type i=0; i<parameters.size(); ++i)
            //  std::cout << parameters[i] << "\t";
            //std::cout << "\n";

            beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.0);
            beagle::real_2d_function_ptr_t rate = drift;

            beagle::real_function_ptr_t localVol
                   = beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(m_Strikes, parameters);
            beagle::real_2d_function_ptr_t volatility
                   = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                   [=](double time, double spot){ return localVol->value(spot); });

            beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                     m_Forward,
                                                                     m_Discounting,
                                                                     drift,
                                                                     volatility,
                                                                     rate,
                                                                     m_Settings);
            auto pODFP = dynamic_cast<beagle::valuation::mixins::OneDimFokkerPlanck*>(forwardPricer.get());

            beagle::dbl_vec_t density(m_Density);
            pODFP->evolve(m_Start, m_End, m_StateVars, density);

            int numStateVars = m_StateVars.size();
            double stateVarStep = (m_StateVars.back() - m_StateVars.front()) / (numStateVars - 1);

            beagle::dbl_vec_t results;
            double fwd = m_Forward->value(m_End);
            for (beagle::dbl_vec_t::size_type j=0; j<m_Strikes.size(); ++j)
            {
              double strike = m_Strikes[j];
              beagle::payoff_ptr_t payoff = strike > fwd ? beagle::product::option::Payoff::call()
                                                         : beagle::product::option::Payoff::put();
              double option(0.0);
              for (int i=0; i<numStateVars; ++i)
              {
                option += stateVarStep * density[i] * payoff->intrinsicValue(fwd * std::exp(m_StateVars[i]), strike);
              }
            
              results.emplace_back(option - m_Targets[j]);
            }

            std::cout << "\nFor this iteration, the erros are: \n";
            for (beagle::dbl_vec_t::size_type i=0; i<results.size(); ++i)
              std::cout << results[i] << "\t";
            std::cout << "\n";

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

            return result;
          }
        private:
          beagle::real_function_ptr_t m_Forward;
          beagle::real_function_ptr_t m_Discounting;
          beagle::valuation::OneDimFiniteDifferenceSettings m_Settings;
          double m_Start;
          double m_End;
          beagle::dbl_vec_t m_StateVars;
          beagle::dbl_vec_t m_Density;
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
        beagle::real_2d_function_ptr_t drift = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.0);
        beagle::real_2d_function_ptr_t rate = drift;

        int centralIndex = volSmiles.back().second.second.size() / 2;
        beagle::real_2d_function_ptr_t volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(
                                                            volSmiles.back().second.second[centralIndex]);

        beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                 forward,
                                                                 discounting,
                                                                 drift,
                                                                 volatility,
                                                                 rate,
                                                                 settings);
        auto pODFP = dynamic_cast<beagle::valuation::mixins::OneDimFokkerPlanck*>(forwardPricer.get());

        beagle::dbl_vec_t stateVars;
        beagle::dbl_vec_t density;
        pODFP->formInitialCondition(volSmiles.back().first, stateVars, density);

        // Now perform calibration by bootstrapping
        beagle::dbl_vec_t expiries;
        beagle::real_function_ptr_coll_t localVolFuncs;
        beagle::dbl_vec_t guesses(volSmiles.front().second.second);

        double start = 0.;
        for (beagle::dbl_vec_t::size_type i=0; i<volSmiles.size(); ++i)
        {
          double end = volSmiles[i].first;
          beagle::dbl_vec_t strikes = volSmiles[i].second.first;

          beagle::calibration_bound_constraint_coll_t constraints(strikes.size(),
                                                                  beagle::calibration::CalibrationBoundConstraint::lowerBoundCalibrationConstraint(0.0001));
          beagle::int_vec_t elimIndices(0U);

          guesses = beagle::calibration::util::getTransformedParameters( guesses, constraints );
          Eigen::VectorXd calibParams(guesses.size());
          for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
            calibParams(i) = guesses[i];

          beagle::calibration_adapter_ptr_t adapter = std::make_shared<LocalVolatilityCalibrationAdapter>(forward,
                                                                                                          discounting,
                                                                                                          settings,
                                                                                                          start,
                                                                                                          end,
                                                                                                          stateVars,
                                                                                                          density,
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
          beagle::real_function_ptr_t localVol
                   = beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(strikes, guesses);
          beagle::real_2d_function_ptr_t volatility
                   = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                   [=](double time, double spot){ return localVol->value(spot); });
          beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEEuroOptionPricer(
                                                                   forward,
                                                                   discounting,
                                                                   drift,
                                                                   volatility,
                                                                   rate,
                                                                   settings);
          auto pODFP = dynamic_cast<beagle::valuation::mixins::OneDimFokkerPlanck*>(forwardPricer.get());
          pODFP->evolve(start, end, stateVars, density);

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