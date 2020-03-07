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
        for (beagle::dbl_vec_t::size_type i=0; i<params.size(); ++i)
          params[i] = constraints[i]->inverseTransform( transformedParamters[i] );

        return params;
      }
      dbl_vec_t getTransformedParameters( const dbl_vec_t& parameters,
                                          const calibration_bound_constraint_coll_t& constraints )
      {
        if (parameters.size() != constraints.size())
          throw("Mismatch in the number of parameters and calibration bound constraints!");

        dbl_vec_t params(parameters);
        for (beagle::dbl_vec_t::size_type i=0; i<params.size(); ++i)
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
      struct PricerAdapter : public CalibrationAdapter
      {
        PricerAdapter( const pricer_ptr_t& pricer,
                       const product_ptr_coll_t& products ) :
          m_Pricer(pricer),
          m_Products(products)
        { }
        virtual ~PricerAdapter( void )
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

          dbl_vec_t result(m_Products.size());
          std::transform(m_Products.cbegin(),
                         m_Products.cend(),
                         result.begin(),
                         [&pricer](const product_ptr_t& product)
                         { return pricer->value(product); });

          return result;
        }
        dbl_mat_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const override
        {
          auto pCWNMP = dynamic_cast<beagle::valuation::mixins::CloneWithNewModelParameters*>(m_Pricer.get());
          if (!pCWNMP)
            throw("Cannot update model parameters");

          dbl_mat_t result(m_Products.size());
          for (auto row : result)
            row.resize(parameters.size());

          double bump = 1e-4;
          for (beagle::dbl_vec_t::size_type j=0; j<parameters.size(); ++j)
          {
            dbl_vec_t forwardParams(parameters);
            dbl_vec_t backwardParams(parameters);
            forwardParams[j] += bump;
            backwardParams[j] -= bump;

            pricer_ptr_t forwardBumpedPricer = pCWNMP->createPricerWithNewModelParameters(forwardParams);
            pricer_ptr_t backwardBumpedPricer = pCWNMP->createPricerWithNewModelParameters(backwardParams);

            for (beagle::dbl_vec_t::size_type i=0; i<m_Products.size(); ++i)
              result[i][j] = ( forwardBumpedPricer->value(m_Products[i])
                             - backwardBumpedPricer->value(m_Products[i]) ) * .5 / bump;
          }

          return result;
        }
      private:
        pricer_ptr_t m_Pricer;
        product_ptr_coll_t m_Products;
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
          beagle::real_function_ptr_t priceFunc = pFD->finiteDifferenceDetails().interpolation()->formFunction(m_Strikes, prices);

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
          dbl_mat_t result(m_InterpStrikes.size());
          for (beagle::dbl_vec_t::size_type i=0; i<m_InterpStrikes.size(); ++i)
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

            for (beagle::dbl_vec_t::size_type i=0; i<m_InterpStrikes.size(); ++i)
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
      for (beagle::dbl_vec_t::size_type j=0; j<constraints.size(); ++j)
      {
        for (beagle::dbl_vec_t::size_type i=0; i<derivs.size(); ++i)
          derivs[i][j] = constraints[j]->transformDerivative(transformedParamters[j], derivs[i][j]);
      }

      return derivs;
    }

    calibration_adapter_ptr_t CalibrationAdapter::pricerAdapter( const pricer_ptr_t& pricer,
                                                                       const product_ptr_coll_t& products )
    {
      return std::make_shared<impl::PricerAdapter>( pricer, products );
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
