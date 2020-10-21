#include "convertible_calibration.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "calibration_adapter.hpp"
#include "calibration_functor.hpp"
#include "util.hpp"
#include "payoff.hpp"
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
        void formModelParameters(const beagle::real_function_ptr_t& smile,
                                 double intensity,
                                 double spot,
                                 double exponent,
                                 beagle::real_2d_function_ptr_t& drift,
                                 beagle::real_2d_function_ptr_t& volatility,
                                 beagle::real_2d_function_ptr_t& rate)
        {
          drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                          [=](double time, double price){ return intensity * std::pow(price / spot, -exponent); } );
          volatility = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(beagle::dbl_vec_t(1U, 1.0),
                                                                                              beagle::real_function_ptr_coll_t(1U, smile));
          rate = drift;
        }

        struct AndersenBuffumCalibrationAdapter : public CalibrationAdapter
        {
          AndersenBuffumCalibrationAdapter(const beagle::real_function_ptr_t& forward,
                                           const beagle::real_function_ptr_t& discounting,
                                           const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                           double start,
                                           double spot,
                                           double exponent,
                                           const beagle::dbl_vec_t& stateVars,
                                           const beagle::dbl_vec_t& density,
                                           const beagle::volatility_smile_credit_spread_coll_t::value_type& quote,
                                           const beagle::interp_builder_ptr_t& interp) :
            m_Forward(forward),
            m_Discounting(discounting),
            m_Settings(settings),
            m_Start(start),
            m_Spot(spot),
            m_Exponent(exponent),
            m_StateVars(stateVars),
            m_Density(density),
            m_Interp(interp)
          {
            m_End = std::get<0>(quote);
            double df = discounting->value(m_End);
            double fwd = forward->value(m_End);

            const auto& smile = std::get<1>(quote);
            m_Strikes = smile.first;
            const auto& vols = smile.second;
            for (beagle::dbl_vec_t::size_type i=0; i<m_Strikes.size(); ++i)
              m_Targets.emplace_back(df * beagle::util::bsCall(m_Strikes[i], fwd, m_End, vols[i]));

            double spread = std::get<2>(quote);
            m_Targets.emplace_back(df * std::exp(-spread * m_End));

            // std::cout << m_Targets[0] << "    " << m_Targets[1] << std::endl << std::endl;
          }
          virtual ~AndersenBuffumCalibrationAdapter( void ) = default;
        public:
          dbl_vec_t values( const dbl_vec_t& parameters ) const override
          {
            beagle::real_2d_function_ptr_t drift;
            beagle::real_2d_function_ptr_t volatility;
            beagle::real_2d_function_ptr_t rate;
            formModelParameters(m_Interp->formFunction(m_Strikes,
                                                       beagle::dbl_vec_t(parameters.cbegin(),
                                                                         parameters.cend() - 1U)),
                                parameters.back(),
                                m_Spot,
                                m_Exponent,
                                drift,
                                volatility,
                                rate);

            beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
                                                                     m_Forward,
                                                                     m_Discounting,
                                                                     drift,
                                                                     volatility,
                                                                     rate,
                                                                     m_Settings);
            auto pODFP = dynamic_cast<beagle::valuation::mixins::OneDimFokkerPlanck*>(forwardPricer.get());

            beagle::dbl_vec_t density(m_Density);
            pODFP->evolve(m_Start, m_End, m_StateVars, density);

            beagle::payoff_ptr_t call = beagle::product::option::Payoff::call();
            beagle::payoff_ptr_t digitalCall = beagle::product::option::Payoff::digitalCall();
            double fwd = m_Forward->value(m_End);
            double df = m_Discounting->value(m_End);

            int numStateVars = m_StateVars.size();
            double stateVarStep = (m_StateVars.back() - m_StateVars.front()) / (numStateVars - 1);

            beagle::dbl_vec_t results;
            for (beagle::dbl_vec_t::size_type j=0; j<m_Strikes.size(); ++j)
            {
              double value(0.0);
              for (int i=0; i<numStateVars; ++i)
              {
                value += stateVarStep * density[i] * call->intrinsicValue(fwd * std::exp(m_StateVars[i]), m_Strikes[j]);
              }

              results.emplace_back(df*value - m_Targets[j]);
            }

            double bond(.0);
            for (int i=0; i<numStateVars; ++i)
            {
              bond += stateVarStep * density[i] * digitalCall->intrinsicValue(1., 0.);
            }

            results.emplace_back(df*bond - m_Targets.back());
            return results;
          }
          dbl_mat_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const override
          {
            beagle::dbl_vec_t::size_type numRows = m_Strikes.size() + 1U;
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
          double m_Spot;
          double m_Exponent;
          beagle::dbl_vec_t m_StateVars;
          beagle::dbl_vec_t m_Density;
          beagle::interp_builder_ptr_t m_Interp;
          beagle::dbl_vec_t m_Strikes;
          beagle::dbl_vec_t m_Targets;
        };
      }

      beagle::andersen_buffum_curve_pair_t
      createCalibratedAndersenBuffumParameters(const beagle::real_function_ptr_t& forward,
                                               const beagle::real_function_ptr_t& discounting,
                                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                               double exponent,
                                               const beagle::dbl_vec_t& expiries,
                                               const beagle::andersen_buffum_param_t& quotes)
      {
        if (expiries.size() != quotes.size())
          throw(std::string("The number of expiries is not the same as the number of calibrations for Andersen-Buffum!"));

        beagle::volatility_smile_credit_spread_coll_t calibQuotes;
        for (beagle::dbl_vec_t::size_type i=0; i<expiries.size(); ++i)
        {
          double expiry = expiries[i];
          double fwd = forward->value(expiry);
          calibQuotes.emplace_back(std::make_tuple(expiry,
                                                   std::make_pair(beagle::dbl_vec_t(1U, fwd),
                                                                  beagle::dbl_vec_t(1U, quotes[i].first)),
                                                   quotes[i].second));
        }

        return createCalibratedAndersenBuffumParameters(forward,
                                                        discounting,
                                                        settings,
                                                        exponent,
                                                        calibQuotes,
                                                        beagle::math::InterpolationBuilder::piecewiseConstantRight());
      }

      beagle::andersen_buffum_curve_pair_t
      createCalibratedAndersenBuffumParameters(const beagle::real_function_ptr_t& forward,
                                               const beagle::real_function_ptr_t& discounting,
                                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                               double exponent,
                                               const beagle::volatility_smile_credit_spread_coll_t& quotes,
                                               const beagle::interp_builder_ptr_t& interp)
      {
        // Set up the spatial finite difference grid and the initial condition
        beagle::real_2d_function_ptr_t drift;
        beagle::real_2d_function_ptr_t volatility;
        beagle::real_2d_function_ptr_t rate;
        double spot = forward->value(0.0);
        formModelParameters(interp->formFunction(std::get<1>(quotes.back()).first,
                                                 std::get<1>(quotes.back()).second),
                            std::get<2>(quotes.back()),
                            spot,
                            exponent,
                            drift,
                            volatility,
                            rate);

        beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
                                                                 forward,
                                                                 discounting,
                                                                 drift,
                                                                 volatility,
                                                                 rate,
                                                                 settings);
        auto pODFP = dynamic_cast<beagle::valuation::mixins::OneDimFokkerPlanck*>(forwardPricer.get());

        beagle::dbl_vec_t stateVars;
        beagle::dbl_vec_t density;
        pODFP->formInitialCondition(std::get<0>(quotes.back()), stateVars, density);

        // Now perform calibration by bootstrapping
        beagle::dbl_vec_t expiries;
        beagle::real_function_ptr_coll_t smiles;
        beagle::dbl_vec_t spreads;
        double start = 0.;
        for (beagle::dbl_vec_t::size_type i=0; i<quotes.size(); ++i)
        {
          const auto& quote = quotes[i];
          double end = std::get<0>(quote);
          const auto& smile = std::get<1>(quote);
          double spread = std::get<2>(quote);

          const beagle::dbl_vec_t& strikes = smile.first;
          const beagle::dbl_vec_t& vols = smile.second;

          beagle::dbl_vec_t::size_type numParams(vols.size() + 1U);
          beagle::dbl_vec_t guesses(vols);
          guesses.emplace_back(spread);

          beagle::calibration_bound_constraint_coll_t constraints(numParams,
                                                                  beagle::calibration::CalibrationBoundConstraint::twoSidedBoundCalibrationConstraint(0., 1.));
          beagle::int_vec_t elimIndices(0U);

          guesses = beagle::calibration::util::getTransformedParameters( guesses, constraints );
          Eigen::VectorXd calibParams(numParams);
          for (beagle::dbl_vec_t::size_type i=0; i<numParams; ++i)
            calibParams(i) = guesses[i];

          beagle::calibration_adapter_ptr_t adapter = std::make_shared<AndersenBuffumCalibrationAdapter>(forward,
                                                                                                         discounting,
                                                                                                         settings,
                                                                                                         start,
                                                                                                         spot,
                                                                                                         exponent,
                                                                                                         stateVars,
                                                                                                         density,
                                                                                                         quote,
                                                                                                         interp);

          beagle::calibration::CalibrationFunctor functor( beagle::dbl_vec_t(numParams, 0.),
                                                           adapter,
                                                           calibParams,
                                                           constraints,
                                                           elimIndices );
          Eigen::LevenbergMarquardt<beagle::calibration::CalibrationFunctor> lm(functor);
          lm.parameters.xtol = 1.0e-6;
          Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(calibParams);

          for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
            guesses[i] = calibParams(i);

          guesses = beagle::calibration::util::getOriginalParameters( guesses, constraints );

          // std::cout << guesses[0] << "    " << guesses[1] << std::endl << std::endl;

          beagle::real_function_ptr_t thisSmile = interp->formFunction(strikes,
                                                                       beagle::dbl_vec_t(guesses.cbegin(),
                                                                                         guesses.cend() - 1U));
          double thisSpread = guesses.back();

          expiries.emplace_back(end);
          smiles.emplace_back(thisSmile);
          spreads.emplace_back(thisSpread);

          // Evolve density
          formModelParameters(thisSmile,
                              thisSpread,
                              spot,
                              exponent,
                              drift,
                              volatility,
                              rate);
          beagle::pricer_ptr_t forwardPricer = beagle::valuation::Pricer::formOneDimForwardPDEArrowDebreuPricer(
                                                                   forward,
                                                                   discounting,
                                                                   drift,
                                                                   volatility,
                                                                   rate,
                                                                   settings);
          auto pODFP = dynamic_cast<beagle::valuation::mixins::OneDimFokkerPlanck*>(forwardPricer.get());
          pODFP->evolve(start, end, stateVars, density);

          start = end;
        }

        beagle::real_2d_function_ptr_t localVol = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(expiries, smiles);
        beagle::real_function_ptr_t intensity = beagle::math::RealFunction::createPiecewiseConstantRightInterpolatedFunction(expiries, spreads);
        beagle::real_2d_function_ptr_t localIntensity = beagle::math::RealTwoDimFunction::createBinaryFunction(
                                                                [=](double time, double price)
                                                                { return intensity->value(time) * std::pow(price / spot, -exponent); } );

        return std::make_pair(localVol, localIntensity);
      }
    }
  }
}