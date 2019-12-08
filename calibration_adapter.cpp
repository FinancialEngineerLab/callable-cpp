#include "calibration_adapter.hpp"

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
            dbl_vec_t backwardParams(parameters);\
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
  }
}
