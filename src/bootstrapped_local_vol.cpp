#include "real_2d_function.hpp"
#include "payoff.hpp"
#include "pricer.hpp"
#include "interpolation_builder.hpp"
#include "real_function.hpp"
#include "calibration_functor.hpp"
#include "unsupported/Eigen/NonLinearOptimization"

#include <iostream>
#include <fstream>

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct BootstrappedLocalVolatilityFunction : public RealTwoDimFunction
      {
        BootstrappedLocalVolatilityFunction( const beagle::dbl_vec_t& expiries,
                                             const beagle::dbl_vec_vec_t& strikesColl,
                                             const beagle::dbl_vec_vec_t& pricesColl,
                                             const beagle::pricer_ptr_t& forwardPricer,
                                             const beagle::payoff_ptr_t& payoff ) :
          m_Expiries( expiries ),
          m_StrikesColl( strikesColl ),
          m_PricesColl( pricesColl ),
          m_ForwardPricer( forwardPricer ),
          m_Payoff( payoff )
        { }
        virtual ~BootstrappedLocalVolatilityFunction( void )
        { }
      public:
        virtual double value( double argX,
                              double argY ) const override
        {
          if (!m_Func)
            doCalibration();

          return m_Func->value( argX, argY );
        }
      private:
        void doCalibration( void ) const
        {
          beagle::dbl_vec_t logStrikes;
          beagle::dbl_vec_t strikes;
          beagle::dbl_vec_t prices;
          beagle::real_function_ptr_coll_t localVols;

          auto pFD = dynamic_cast<beagle::valuation::mixins::FiniteDifference*>(m_ForwardPricer.get());
          double spot = pFD->spot();
          double rate = pFD->rate();
          int stepsPerAnnum = pFD->stepsPerAnnum();
          int numStrikes = pFD->numberOfStateVariableSteps();
          double numStdDev = pFD->numberOfStandardDeviations();
          const beagle::discrete_dividend_schedule_t& dividends = pFD->dividends();
          const beagle::dividend_policy_ptr_t& policy = pFD->dividendPolicy();
          const beagle::interp_builder_ptr_t& interp = pFD->interpolation();
          pFD->formStateVariableSteps( m_Expiries.back(), logStrikes, strikes );

          auto pOVCP = dynamic_cast<beagle::valuation::mixins::OptionValueCollectionProvider*>(m_ForwardPricer.get());
          if (!pOVCP)
            throw(std::string("Cast for OptionValueCollectionProvider failed!"));
          pOVCP->formInitialOptionValueCollection( m_Payoff, strikes, prices );

          double start(0.0);
          for (beagle::dbl_vec_t::size_type i=0; i<m_Expiries.size(); ++i)
          {
            double end = m_Expiries[i];
            calibration_adapter_ptr_t adapter = beagle::calibration::CalibrationAdapter::forwardPDEPricerAdapter(m_ForwardPricer,
                                                                                                                 start,
                                                                                                                 end,
                                                                                                                 m_Payoff,
                                                                                                                 logStrikes,
                                                                                                                 strikes,
                                                                                                                 prices,
                                                                                                                 m_StrikesColl[i]);

            beagle::dbl_vec_t guesses{.3, .3, .3, .3, .3, .3, .3, .3, .3};
            beagle::calibration_bound_constraint_coll_t constraints(guesses.size(),
                                                                    beagle::calibration::CalibrationBoundConstraint::twoSidedBoundCalibrationConstraint(0., 2.));
            beagle::int_vec_t elimIndices(0U);

            guesses = beagle::calibration::util::getTransformedParameters( guesses, constraints );
            Eigen::VectorXd calibParams(guesses.size());
            for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
              calibParams(i) = guesses[i];

            beagle::calibration::CalibrationFunctor functor( m_PricesColl[i], adapter, calibParams, constraints, elimIndices );
            Eigen::LevenbergMarquardt<beagle::calibration::CalibrationFunctor> lm(functor);
            Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(calibParams);

            for (beagle::dbl_vec_t::size_type i=0; i<guesses.size(); ++i)
              guesses[i] = calibParams(i);

            guesses = beagle::calibration::util::getOriginalParameters( guesses, constraints );

            start = end;
            prices = m_PricesColl[i];
            localVols.push_back( beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(m_StrikesColl[i],
                                                                                                                   guesses) );
          }

          m_Func = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( m_Expiries, localVols );
        }
      private:
        beagle::dbl_vec_t m_Expiries;
        beagle::dbl_vec_vec_t m_StrikesColl;
        beagle::dbl_vec_vec_t m_PricesColl;
        beagle::pricer_ptr_t m_ForwardPricer;
        beagle::payoff_ptr_t m_Payoff;
        mutable beagle::real_2d_function_ptr_t m_Func;
      };

    }

    beagle::real_2d_function_ptr_t
    RealTwoDimFunction::createBootstrappedLocalVolatilityFunction(
                                  const beagle::dbl_vec_t& expiries,
                                  const beagle::dbl_vec_vec_t& strikesColl,
                                  const beagle::dbl_vec_vec_t& pricesColl,
                                  const beagle::pricer_ptr_t& forwardPricer,
                                  const beagle::payoff_ptr_t& payoff )
    {
      if (expiries.size() != strikesColl.size())
        throw(std::string("Number of expiries must be identical to the number of strike collection"));
      if (expiries.size() != pricesColl.size())
        throw(std::string("Number of expiries must be identical to the number of price collection"));
      for (beagle::dbl_vec_vec_t::size_type i=0; i<strikesColl.size(); ++i)
      {
        if (strikesColl[i].size() != pricesColl[i].size())
          throw(std::string("The number of strikes must be identical to the number of prices"));
      }

      return std::make_shared<impl::BootstrappedLocalVolatilityFunction>( expiries,
                                                                          strikesColl,
                                                                          pricesColl,
                                                                          forwardPricer,
                                                                          payoff );
    }
  }
}
