#include "calibration_adapter.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "interpolation_builder.hpp"

namespace beagle
{
  namespace calibration
  {
    namespace util
    {
      dbl_vec_t getOriginalParameters( const dbl_vec_t& transformedParamters,
                                       const calibration_bound_constraint_coll_t& constraints )
      {
        if (transformedParamters.size() != constraints.size())
          throw("Mismatch in the number of transformed parameters and calibration bound constraints!");

        dbl_vec_t params(transformedParamters);
        for (int i=0; i<params.size(); ++i)
          params[i] = constraints[i]->inverseTransform( transformedParamters[i] );

        return params;
      }
      dbl_vec_t getTransformedParameters( const dbl_vec_t& parameters,
                                          const calibration_bound_constraint_coll_t& constraints )
      {
        if (parameters.size() != constraints.size())
          throw("Mismatch in the number of parameters and calibration bound constraints!");

        dbl_vec_t params(parameters);
        for (int i=0; i<params.size(); ++i)
          params[i] = constraints[i]->transform( parameters[i] );

        return params;
      }
      dbl_vec_t liveValues( const calibration_adapter_ptr_t& adapter,
                            const dbl_vec_t& parameters )
      {
        return adapter->values(parameters);
      }
      dbl_vec_t liveValues( const calibration_adapter_ptr_t& adapter,
                            const dbl_vec_t& transformedParamters,
                            const calibration_bound_constraint_coll_t& constraints )
      {
        return adapter->values(util::getOriginalParameters(transformedParamters, constraints));
      }
      dbl_mat_t jacobianMatrix( const calibration_adapter_ptr_t& adapter,
                                const dbl_vec_t& parameters )
      {
        return adapter->derivativeWithRespectToParameters( parameters );
      }

      dbl_mat_t jacobianMatrix( const calibration_adapter_ptr_t& adapter,
                                const dbl_vec_t& transformedParamters,
                                const calibration_bound_constraint_coll_t& constraints )
      {
        return adapter->derivativeWithRespectToParameters( util::getOriginalParameters(transformedParamters, constraints) );
      }
    }

    namespace impl
    {
      struct OptionPricerAdapter : public CalibrationAdapter
      {
        OptionPricerAdapter( const pricer_ptr_t& pricer,
                             const option_ptr_coll_t& options ) :
          m_Pricer(pricer),
          m_Options(options)
        { }
        virtual ~OptionPricerAdapter( void )
        { }
      public:
        dbl_vec_t values( const dbl_vec_t& parameters ) const override
        {
          pricer_ptr_t pricer;
          auto pCWNMP = dynamic_cast<beagle::valuation::mixins::CloneWithNewModelParameters*>(m_Pricer.get());
          if (pCWNMP)
            pricer = pCWNMP->createPricerWithNewModelParameters(parameters);
          else
            throw("Cannot update model parameters");

          dbl_vec_t result(m_Options.size());
          std::transform(m_Options.cbegin(),
                         m_Options.cend(),
                         result.begin(),
                         [&pricer](const option_ptr_t& option)
                         { return pricer->optionValue(option); });

          return result;
        }
        dbl_mat_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const override
        {
          auto pCWNMP = dynamic_cast<beagle::valuation::mixins::CloneWithNewModelParameters*>(m_Pricer.get());
          if (!pCWNMP)
            throw("Cannot update model parameters");

          dbl_mat_t result(m_Options.size());
          for (auto row : result)
            row.resize(parameters.size());

          double bump = 1e-4;
          for (int j=0; j<parameters.size(); ++j)
          {
            dbl_vec_t forwardParams(parameters);
            dbl_vec_t backwardParams(parameters);
            forwardParams[j] += bump;
            backwardParams[j] -= bump;

            pricer_ptr_t forwardBumpedPricer = pCWNMP->createPricerWithNewModelParameters(forwardParams);
            pricer_ptr_t backwardBumpedPricer = pCWNMP->createPricerWithNewModelParameters(backwardParams);

            for (int i=0; i<m_Options.size(); ++i)
              result[i][j] = ( forwardBumpedPricer->optionValue(m_Options[i])
                             - backwardBumpedPricer->optionValue(m_Options[i]) ) * .5 / bump;
          }

          return result;
        }
      private:
        pricer_ptr_t m_Pricer;
        option_ptr_coll_t m_Options;
      };

      struct ForwardPDEPricerAdapter : public CalibrationAdapter
      {
        ForwardPDEPricerAdapter( const pricer_ptr_t& forwardPricer,
                                 double start,
                                 double end,
                                 const beagle::payoff_ptr_t& payoff,
                                 const beagle::dbl_vec_t& logStrikes,
                                 const beagle::dbl_vec_t& strikes,
                                 const beagle::dbl_vec_t& prices,
                                 const beagle::dbl_vec_t& interpStrikes ) :
          m_Pricer(forwardPricer),
          m_Start(start),
          m_End(end),
          m_Payoff(payoff),
          m_LogStrikes(logStrikes),
          m_Strikes(strikes),
          m_Prices(prices),
          m_InterpStrikes(interpStrikes)
        { }
        virtual ~ForwardPDEPricerAdapter( void )
        { }
      public:
        dbl_vec_t values( const dbl_vec_t& parameters ) const override
        {
          pricer_ptr_t pricer;
          auto pCWNMP = dynamic_cast<beagle::valuation::mixins::CloneWithNewLocalVolatilitySurface*>(m_Pricer.get());
          if (pCWNMP)
          {
            beagle::real_2d_function_ptr_t newLocalVolSurface
              = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(
                        dbl_vec_t(1U, 1.0),
                        real_function_ptr_coll_t(
                              1U,
                              beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(m_InterpStrikes,
                                                                                                                parameters)));
            pricer = pCWNMP->createPricerWithNewLocalVolatilitySurface(newLocalVolSurface);
          }
          else
          {
            throw("Cannot update model parameters");
          }

          dbl_vec_t prices(m_Prices);
          auto pOVCP = dynamic_cast<beagle::valuation::mixins::OptionValueCollectionProvider*>(pricer.get());
          pOVCP->optionValueCollection(m_Start,
                                       m_End,
                                       m_Payoff,
                                       m_LogStrikes,
                                       m_Strikes,
                                       prices );

          auto pFD = dynamic_cast<beagle::valuation::mixins::FiniteDifference*>(pricer.get());
          beagle::real_function_ptr_t priceFunc = pFD->interpolation()->formFunction(m_Strikes, prices);

          dbl_vec_t result(m_InterpStrikes.size());
          std::transform(m_InterpStrikes.cbegin(),
                         m_InterpStrikes.cend(),
                         result.begin(),
                         [&priceFunc](double strike)
                         { return priceFunc->value(strike); });
          return result;
        }
        dbl_mat_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const override
        {
          auto pCWNMP = dynamic_cast<beagle::valuation::mixins::CloneWithNewLocalVolatilitySurface*>(m_Pricer.get());
          if (!pCWNMP)
            throw("Cannot update model parameters");

          dbl_mat_t result(m_InterpStrikes.size());
          for (int i=0; i<m_InterpStrikes.size(); ++i)
          {
            result[i].resize(parameters.size());
          }

          double bump = 1e-4;
          for (int j=0; j<parameters.size(); ++j)
          {
            dbl_vec_t forwardParams(parameters);
            dbl_vec_t backwardParams(parameters);
            forwardParams[j] += bump;
            backwardParams[j] -= bump;

            beagle::real_2d_function_ptr_t forwardSurface
              = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(
                        dbl_vec_t(1U, 1.0),
                        real_function_ptr_coll_t(
                              1U,
                              beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(m_InterpStrikes,
                                                                                                                forwardParams)));
            beagle::real_2d_function_ptr_t backwardSurface
              = beagle::math::RealTwoDimFunction::createPiecewiseConstantRightFunction(
                        dbl_vec_t(1U, 1.0),
                        real_function_ptr_coll_t(
                              1U,
                              beagle::math::RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(m_InterpStrikes,
                                                                                                                backwardParams)));

            pricer_ptr_t forwardBumpedPricer = pCWNMP->createPricerWithNewLocalVolatilitySurface(forwardSurface);
            pricer_ptr_t backwardBumpedPricer = pCWNMP->createPricerWithNewLocalVolatilitySurface(backwardSurface);

            dbl_vec_t forwardPrices(m_Prices);
            auto pOVCP = dynamic_cast<beagle::valuation::mixins::OptionValueCollectionProvider*>(forwardBumpedPricer.get());
            pOVCP->optionValueCollection(m_Start,
                                         m_End,
                                         m_Payoff,
                                         m_LogStrikes,
                                         m_Strikes,
                                         forwardPrices);

            auto pFD = dynamic_cast<beagle::valuation::mixins::FiniteDifference*>(forwardBumpedPricer.get());
            beagle::real_function_ptr_t priceFunc = pFD->interpolation()->formFunction(m_Strikes, forwardPrices);

            dbl_vec_t forwardResult(m_InterpStrikes.size());
            std::transform(m_InterpStrikes.cbegin(),
                           m_InterpStrikes.cend(),
                           forwardResult.begin(),
                           [&priceFunc](double strike)
                           { return priceFunc->value(strike); });

            dbl_vec_t backwardPrices(m_Prices);
            auto pOVCQ = dynamic_cast<beagle::valuation::mixins::OptionValueCollectionProvider*>(backwardBumpedPricer.get());
            pOVCQ->optionValueCollection(m_Start,
                                         m_End,
                                         m_Payoff,
                                         m_LogStrikes,
                                         m_Strikes,
                                         backwardPrices);

            auto pFE = dynamic_cast<beagle::valuation::mixins::FiniteDifference*>(backwardBumpedPricer.get());
            beagle::real_function_ptr_t priceFund = pFE->interpolation()->formFunction(m_Strikes, backwardPrices);

            dbl_vec_t backwardResult(m_InterpStrikes.size());
            std::transform(m_InterpStrikes.cbegin(),
                           m_InterpStrikes.cend(),
                           backwardResult.begin(),
                           [&priceFund](double strike)
                           { return priceFund->value(strike); });

            for (int i=0; i<m_InterpStrikes.size(); ++i)
            {
              result[i][j] = ( forwardResult[i] - backwardResult[i] ) * .5 / bump;
            }
          }

          return result;
        }
      private:
        pricer_ptr_t m_Pricer;
        double m_Start;
        double m_End;
        beagle::payoff_ptr_t m_Payoff;
        beagle::dbl_vec_t m_LogStrikes;
        beagle::dbl_vec_t m_Strikes;
        beagle::dbl_vec_t m_Prices;
        beagle::dbl_vec_t m_InterpStrikes;
      };
    }

    CalibrationAdapter::~CalibrationAdapter( void )
    { }

    dbl_vec_t CalibrationAdapter::values( const dbl_vec_t& transformedParamters,
                                          const calibration_bound_constraint_coll_t& constraints ) const
    {
      return values( util::getOriginalParameters(transformedParamters, constraints) );
    }

    dbl_mat_t CalibrationAdapter::derivativeWithRespectToTransformedParameters(
                     const dbl_vec_t& transformedParamters,
                     const calibration_bound_constraint_coll_t& constraints ) const
    {
      dbl_mat_t derivs = derivativeWithRespectToParameters( util::getOriginalParameters(transformedParamters, constraints) );
      for (int j=0; j<constraints.size(); ++j)
      {
        for (int i=0; i<derivs.size(); ++i)
          derivs[i][j] = constraints[j]->transformDerivative(transformedParamters[j], derivs[i][j]);
      }

      return derivs;
    }

    calibration_adapter_ptr_t CalibrationAdapter::optionPricerAdapter( const pricer_ptr_t& pricer,
                                                                       const option_ptr_coll_t& options )
    {
      return std::make_shared<impl::OptionPricerAdapter>( pricer, options );
    }

    calibration_adapter_ptr_t CalibrationAdapter::forwardPDEPricerAdapter( const pricer_ptr_t& forwardPricer,
                                                                           double start,
                                                                           double end,
                                                                           const beagle::payoff_ptr_t& payoff,
                                                                           const beagle::dbl_vec_t& logStrikes,
                                                                           const beagle::dbl_vec_t& strikes,
                                                                           const beagle::dbl_vec_t& prices,
                                                                           const beagle::dbl_vec_t& interpStrikes )
    {
      return std::make_shared<impl::ForwardPDEPricerAdapter>( forwardPricer,
                                                              start,
                                                              end,
                                                              payoff,
                                                              logStrikes,
                                                              strikes,
                                                              prices,
                                                              interpStrikes );
    }
  }
}
