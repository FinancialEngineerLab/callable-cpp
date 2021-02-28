#include "closed_form_calibration.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "calibration_adapter.hpp"
#include "calibration_functor.hpp"
#include "util.hpp"
#include "payoff.hpp"
#include "option.hpp"
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
        struct ClosedFormEuropeanOptionPricersAdapter : public CalibrationAdapter
        {
          ClosedFormEuropeanOptionPricersAdapter(const beagle::pricer_ptr_t& pricer,
                                                 const beagle::product_ptr_coll_t& products) :
            m_Pricer(pricer),
            m_Products(products)
          { }
        public:
          beagle::dbl_vec_t values( const beagle::dbl_vec_t& parameters ) const override
          {
            beagle::real_function_ptr_coll_t funcs;
            funcs.reserve(parameters.size());
            for (double p : parameters)
              funcs.push_back(beagle::math::RealFunction::createConstantFunction(p));

            auto pCF = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(m_Pricer.get());
            beagle::pricer_ptr_t pricer = pCF->updateModelParameters(funcs);

            beagle::dbl_vec_t results;
            results.reserve(m_Products.size());
            for (const auto& p : m_Products)
              results.push_back(pricer->value(p));

            return results;
          }
          beagle::dbl_mat_t derivativeWithRespectToParameters( const beagle::dbl_vec_t& parameters ) const override
          {
            int numProducts = m_Products.size();
            int numParams = parameters.size();
            beagle::dbl_mat_t jacobian(numProducts, beagle::dbl_vec_t(numParams, 0.));

            beagle::dbl_vec_t bumpUpParams(parameters);
            beagle::dbl_vec_t bumpDownParams(parameters);
            double bumpSize = 1e-5;
            double invBumpSize = 1. / bumpSize;
            for (int i = 0; i < numParams; ++i)
            {
              bumpUpParams[i] += bumpSize;
              bumpDownParams[i] -= bumpSize;

              beagle::dbl_vec_t bumpUpValues = values(bumpUpParams);
              beagle::dbl_vec_t bumpDownValues = values(bumpDownParams);

              for (int j = 0; j < numProducts; ++j)
                jacobian[j][i] = (bumpUpValues[j] - bumpDownValues[j]) * .5 * invBumpSize;

              bumpUpParams[i] = parameters[i];
              bumpDownParams[i] = parameters[i];
            }

            return jacobian;
          }
        private:
          beagle::pricer_ptr_t m_Pricer;
          beagle::product_ptr_coll_t m_Products;
        };
      }

      beagle::real_function_ptr_coll_t
      createCalibratedClosedFormEuropeanOptionPricerParameters(const beagle::real_function_ptr_t& forward,
                                                               const beagle::real_function_ptr_t& discounting,
                                                               const beagle::volatility_smile_coll_t& volSmiles,
                                                               const beagle::pricer_ptr_t& pricer,
                                                               const beagle::dbl_vec_t& guesses,
                                                               const beagle::calibration_bound_constraint_coll_t& constraints,
                                                               const beagle::interp_builder_ptr_coll_t& interps,
                                                               double tolerance)
      {
        auto pCF = dynamic_cast<beagle::valuation::mixins::ClosedFormEuroOption*>(pricer.get());
        int numParams = pCF->numberOfParameters();
        if (guesses.size() != numParams)
          throw(std::string("The number of initial guesses must agree with the number of calibration parameters"));
        if (constraints.size() != guesses.size())
          throw(std::string("The number of calibration bound constraints is different from the number of initial guesses"));
        if (interps.size() != guesses.size())
          throw(std::string("The number of interpolation methods is different from the number of initial guesses"));

        beagle::payoff_ptr_t payoff = beagle::product::option::Payoff::call();

        int numExpiries = volSmiles.size();
        beagle::dbl_mat_t results(numParams, beagle::dbl_vec_t(numExpiries, 0.0));

        beagle::dbl_vec_t expiries(numExpiries, 0.);
        for (int i = 0; i < numExpiries; ++i)
        {
          double expiry = volSmiles[i].first;
          expiries[i] = expiry;

          const beagle::dbl_vec_t& strikes = volSmiles[i].second.first;
          const beagle::dbl_vec_t& vols = volSmiles[i].second.second;

          beagle::product_ptr_coll_t options(strikes.size());
          beagle::dbl_vec_t targets(strikes.size());
          for (int j = 0; j < strikes.size(); ++j)
          {
            options[j] = beagle::product::option::Option::createEuropeanOption(expiry, strikes[j], payoff);
            targets[j] = beagle::valuation::Pricer::formBlackScholesClosedFormEuropeanOptionPricer(forward, discounting, vols[j])->value(options[j]);
          }

          beagle::calibration_adapter_ptr_t adapter = std::make_shared<ClosedFormEuropeanOptionPricersAdapter>(pricer, options);

          beagle::dbl_vec_t calibResults = beagle::calibration::util::getTransformedParameters( guesses, constraints );
          Eigen::VectorXd calibParams(numParams);
          for (int j=0; j < numParams; ++j)
            calibParams(j) = calibResults[j];

          beagle::calibration::CalibrationFunctor functor(targets,
                                                          adapter,
                                                          calibParams,
                                                          constraints,
                                                          {} );

          Eigen::LevenbergMarquardt<beagle::calibration::CalibrationFunctor> lm(functor);
          lm.parameters.xtol = tolerance;
          Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(calibParams);

          for (int j=0; j < numParams; ++j)
            calibResults[j] = calibParams(j);

          calibResults = beagle::calibration::util::getOriginalParameters( calibResults, constraints );
          for (int j = 0; j < numParams; ++j)
            results[j][i] = calibResults[j];
        }

        beagle::real_function_ptr_coll_t calibFuncs(numParams);
        for (int i = 0; i < numParams; ++i)
          calibFuncs[i] = interps[i]->formFunction(expiries, results[i]);

        return calibFuncs;
      }
    }
  }
}