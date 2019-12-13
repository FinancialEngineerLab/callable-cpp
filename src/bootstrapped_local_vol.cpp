#include "real_2d_function.hpp"
#include "payoff.hpp"
#include "pricer.hpp"
#include "interpolation_builder.hpp"
#include "real_function.hpp"
#include "calibration_functor.hpp"
#include "unsupported/Eigen/NonLinearOptimization"

#include <iostream>

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct BootstrappedLocalVolatilityFunction : public RealTwoDimFunction
      {
        BootstrappedLocalVolatilityFunction( const beagle::dbl_vec_t& expiries,
                                             const beagle::dbl_vec_t& initialGuesses,
                                             const beagle::dbl_vec_vec_t& strikesColl,
                                             const beagle::dbl_vec_vec_t& pricesColl,
                                             const beagle::pricer_ptr_t& forwardPricer,
                                             const beagle::payoff_ptr_t& payoff,
                                             const beagle::interp_builder_ptr_t& interp ) :
          m_Expiries( expiries ),
          m_InitialGuesses( initialGuesses ),
          m_StrikesColl( strikesColl ),
          m_PricesColl( pricesColl ),
          m_ForwardPricer( forwardPricer ),
          m_Payoff( payoff ),
          m_Interpolation( interp )
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
          const beagle::valuation::FiniteDifferenceDetails& fdDetails = pFD->finiteDifferenceDetails();
          fdDetails.formStateVariableSteps( m_Expiries.back(), logStrikes, strikes );

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

            beagle::dbl_vec_t guesses(m_StrikesColl[i].size(), m_InitialGuesses[i]);
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
            beagle::real_function_ptr_t localVol = m_Interpolation->formFunction(m_StrikesColl[i], guesses);
            beagle::pricer_ptr_t forwardPricer
              = beagle::valuation::Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer(
                                                fdDetails,
                                                beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(
                                                        dbl_vec_t(1U, 1.0),
                                                        beagle::real_function_ptr_coll_t(1U, localVol)));

            auto pOVCP = dynamic_cast<beagle::valuation::mixins::OptionValueCollectionProvider*>(forwardPricer.get());
            pOVCP->optionValueCollection(start, end, m_Payoff, logStrikes, strikes, prices);

            start = end;
            localVols.push_back( localVol );
          }

          m_Func = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction( m_Expiries, localVols );
        }
      private:
        beagle::dbl_vec_t m_Expiries;
        beagle::dbl_vec_t m_InitialGuesses;
        beagle::dbl_vec_vec_t m_StrikesColl;
        beagle::dbl_vec_vec_t m_PricesColl;
        beagle::pricer_ptr_t m_ForwardPricer;
        beagle::payoff_ptr_t m_Payoff;
        beagle::interp_builder_ptr_t m_Interpolation;
        mutable beagle::real_2d_function_ptr_t m_Func;
      };

    }

    beagle::real_2d_function_ptr_t
    RealTwoDimFunction::createBootstrappedLocalVolatilityFunction(
                                  const beagle::dbl_vec_t& expiries,
                                  const beagle::dbl_vec_t& initialGuesses,
                                  const beagle::dbl_vec_vec_t& strikesColl,
                                  const beagle::dbl_vec_vec_t& pricesColl,
                                  const beagle::pricer_ptr_t& forwardPricer,
                                  const beagle::payoff_ptr_t& payoff,
                                  const beagle::interp_builder_ptr_t& interp )
    {
      if (expiries.size() != initialGuesses.size())
        throw(std::string("Number of expiries must be identical to the number of initial guesses"));
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
                                                                          initialGuesses,
                                                                          strikesColl,
                                                                          pricesColl,
                                                                          forwardPricer,
                                                                          payoff,
                                                                          interp );
    }
  }
}
