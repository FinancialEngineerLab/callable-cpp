#include "convertible_calibration.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "calibration_adapter.hpp"
#include "calibration_functor.hpp"
#include "util.hpp"
#include "payoff.hpp"
#include "unsupported/Eigen/NonLinearOptimization"

namespace beagle
{
  namespace calibration
  {
    namespace util
    {
      namespace
      {
        void formModelParameters(double sigma,
                                 double intensity,
                                 double spot,
                                 double exponent,
                                 beagle::real_2d_function_ptr_t& drift,
                                 beagle::real_2d_function_ptr_t& volatility,
                                 beagle::real_2d_function_ptr_t& rate)
        {
          drift = beagle::math::RealTwoDimFunction::createBinaryFunction(
                          [=](double time, double price){ return intensity * std::pow(price / spot, -exponent); } );
          volatility = beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(sigma);
          rate = drift;
        }
        
        struct AndersenBuffumCalibrationAdapter : public CalibrationAdapter
        {
          AndersenBuffumCalibrationAdapter(const beagle::real_function_ptr_t& forward,
                                           const beagle::real_function_ptr_t& discounting,
                                           const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                           double start,
                                           double end,
                                           double spot,
                                           double exponent,
                                           const beagle::dbl_vec_t& stateVars,
                                           const beagle::dbl_vec_t& density,
                                           const beagle::two_dbl_t& quotes) :
            m_Forward(forward),
            m_Discounting(discounting),
            m_Settings(settings),
            m_Start(start),
            m_End(end),
            m_Spot(spot),
            m_Exponent(exponent),
            m_StateVars(stateVars),
            m_Density(density)
          {
            double df = discounting->value(end);
            double strike = forward->value(end);
            m_Targets.first = df * beagle::util::bsCall(strike, strike, end, quotes.first);
            m_Targets.second = df * std::exp(-quotes.second * end);
          }
          virtual ~AndersenBuffumCalibrationAdapter( void )
          { }
        public:
          dbl_vec_t values( const dbl_vec_t& parameters ) const override
          {
            if (parameters.size() != 2U)
              throw(std::string("Can only calibration exactly two parameters!"));

            beagle::real_2d_function_ptr_t drift;
            beagle::real_2d_function_ptr_t volatility;
            beagle::real_2d_function_ptr_t rate;
            formModelParameters(parameters[0],
                                parameters[1],
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
            double strike = m_Forward->value(m_End);
            double df = m_Discounting->value(m_End);

            int numStateVars = m_StateVars.size();
            double stateVarStep = (m_StateVars.back() - m_StateVars.front()) / (numStateVars - 1);

            double option(.0);
            double bond(.0);
            for (int i=0; i<numStateVars; ++i)
            {
              option += stateVarStep * density[i] * call->intrinsicValue(strike * std::exp(m_StateVars[i]), strike);
              bond += stateVarStep * density[i] * digitalCall->intrinsicValue(strike * std::exp(m_StateVars[i]), 0.);
            }

            return beagle::dbl_vec_t{df * option - m_Targets.first,
                                     df * bond - m_Targets.second};
          }
          dbl_mat_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const override
          {
            beagle::dbl_vec_t::size_type numRows = 2U;
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
          beagle::two_dbl_t m_Targets;
        };
      }

      beagle::andersen_buffum_param_t
      createCalibratedAndersenBuffumParameters(const beagle::real_function_ptr_t& forward,
                                               const beagle::real_function_ptr_t& discounting,
                                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings,
                                               double exponent,
                                               const beagle::dbl_vec_t& expiries,
                                               const beagle::andersen_buffum_param_t& quotes)
      {
        if (expiries.size() != quotes.size())
          throw(std::string("The number of expiries is not the same as the number of calibrations for Andersen-Buffum!"));

        // Set up the spatial finite difference grid and the initial condition
        beagle::real_2d_function_ptr_t drift;
        beagle::real_2d_function_ptr_t volatility;
        beagle::real_2d_function_ptr_t rate;
        double spot = forward->value(0.0);
        formModelParameters(quotes.back().first,
                            quotes.back().second,
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
        pODFP->formInitialCondition(expiries.back(), stateVars, density);

        beagle::dbl_vec_t::size_type numParams(2U);

        // Now perform calibration by bootstrapping
        beagle::andersen_buffum_param_t params;
        double start = 0.;
        for (beagle::dbl_vec_t::size_type i=0; i<expiries.size(); ++i)
        {
          double end = expiries[i];
          beagle::dbl_vec_t guesses{quotes[i].first, quotes[i].second};

          beagle::calibration_bound_constraint_coll_t constraints(numParams,
                                                                  beagle::calibration::CalibrationBoundConstraint::lowerBoundCalibrationConstraint(0.));
          beagle::int_vec_t elimIndices(0U);

          guesses = beagle::calibration::util::getTransformedParameters( guesses, constraints );
          Eigen::VectorXd calibParams(guesses.size());
          for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
            calibParams(i) = guesses[i];

          beagle::calibration_adapter_ptr_t adapter = std::make_shared<AndersenBuffumCalibrationAdapter>(forward,
                                                                                                         discounting,
                                                                                                         settings,
                                                                                                         start,
                                                                                                         end,
                                                                                                         spot,
                                                                                                         exponent,
                                                                                                         stateVars,
                                                                                                         density,
                                                                                                         quotes[i]);

          beagle::calibration::CalibrationFunctor functor( beagle::dbl_vec_t(2U, 0.),
                                                           adapter,
                                                           calibParams,
                                                           constraints,
                                                           elimIndices );
          Eigen::LevenbergMarquardt<beagle::calibration::CalibrationFunctor> lm(functor);
          lm.parameters.xtol = 1.0e-4;
          Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(calibParams);

          for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
            guesses[i] = calibParams(i);

          guesses = beagle::calibration::util::getOriginalParameters( guesses, constraints );
          params.emplace_back(guesses[0], guesses[1]);

          // Evolve density
          formModelParameters(guesses[0],
                              guesses[1],
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

          start = expiries[i];
        }

        return params;
      }
    }
  }
}