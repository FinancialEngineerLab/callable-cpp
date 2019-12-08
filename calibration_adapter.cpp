#include "calibration_adapter.hpp"

namespace sabr_test
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
          params(i) = constraints[i]->inverseTransform( transformedParamters(i) );

        return params;
      }
      dbl_vec_t getTransformedParameters( const dbl_vec_t& parameters,
                                          const calibration_bound_constraint_coll_t& constraints )
      {
        if (parameters.size() != constraints.size())
          throw("Mismatch in the number of parameters and calibration bound constraints!");

        dbl_vec_t params(parameters);
        for (int i=0; i<params.size(); ++i)
          params(i) = constraints[i]->transform( parameters(i) );

        return params;
      }
      dbl_vec_t liveValues( const calibration_adapter_coll_t& adapters,
                            const dbl_vec_t& parameters )
      {
        dbl_vec_t liveValues( adapters.size() );
        for (int i=0; i<liveValues.size(); ++i)
          liveValues(i) = adapters[i]->value(parameters);

        return liveValues;
      }
      dbl_vec_t liveValues( const calibration_adapter_coll_t& adapters,
                            const dbl_vec_t& transformedParamters,
                            const calibration_bound_constraint_coll_t& constraints )
      {
        if (transformedParamters.size() != constraints.size())
          throw("Mismatch in the number of transformed parameters and calibration bound constraints!");

        dbl_vec_t liveValues( adapters.size() );
        for (int i=0; i<liveValues.size(); ++i)
          liveValues(i) = adapters[i]->value(transformedParamters, constraints);

        return liveValues;
      }
      dbl_mat_t jacobianMatrix( const calibration_adapter_coll_t& adapters,
                                const dbl_vec_t& parameters )
      {
        dbl_mat_t jacobian( adapters.size(), parameters.size() );
        for (int i=0; i<jacobian.rows(); ++i)
        {
          for (int j=0; j<jacobian.cols(); ++j)
          {
            dbl_vec_t derivs = adapters[i]->derivativeWithRespectToParameters( parameters );
            jacobian(i, j) = derivs(j);
          }
        }

        return jacobian;
      }

      dbl_mat_t jacobianMatrix( const calibration_adapter_coll_t& adapters,
                                const dbl_vec_t& transformedParamters,
                                const calibration_bound_constraint_coll_t& constraints )
      {
        if (transformedParamters.size() != constraints.size())
          throw("Mismatch in the number of transformed parameters and calibration bound constraints!");

        dbl_mat_t jacobian( adapters.size(), transformedParamters.size() );
        for (int i=0; i<jacobian.rows(); ++i)
        {
          for (int j=0; j<jacobian.cols(); ++j)
          {
            dbl_vec_t derivs = adapters[i]->derivativeWithRespectToTransformedParameters( transformedParamters, constraints );
            jacobian(i, j) = derivs(j);
          }
        }

        return jacobian;
      }
    }

    namespace impl
    {
      struct OptionPricerAdapter : public CalibrationAdapter
      {
        OptionPricerAdapter( const opt_ptr_t& optionPricer,
                             double strike,
                             double forward,
                             double expiry ) :
          m_Option(optionPricer),
          m_Strike(strike),
          m_Forward(forward),
          m_Expiry(expiry)
        { }
        virtual ~OptionPricerAdapter( void )
        { }
      public:
        dbl_t value( const dbl_vec_t& parameters ) const override
        {
          return m_Option->updateModelParameters( parameters )->callValue(m_Strike, m_Forward, m_Expiry);
        }
        dbl_vec_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const override
        {
          return m_Option->updateModelParameters( parameters )->derivativeWithRespectToParameters(m_Strike, m_Forward, m_Expiry);
        }
      private:
        opt_ptr_t m_Option;
        dbl_t m_Strike;
        dbl_t m_Forward;
        dbl_t m_Expiry;
      };
    }

    CalibrationAdapter::~CalibrationAdapter( void )
    { }

    dbl_t CalibrationAdapter::value( const dbl_vec_t& transformedParamters,
                                     const calibration_bound_constraint_coll_t& constraints ) const
    {
      return value( util::getOriginalParameters(transformedParamters, constraints) );
    }

    dbl_vec_t CalibrationAdapter::derivativeWithRespectToTransformedParameters(
                     const dbl_vec_t& transformedParamters,
                     const calibration_bound_constraint_coll_t& constraints ) const
    {
      dbl_vec_t derivs = derivativeWithRespectToParameters( util::getOriginalParameters(transformedParamters, constraints) );
      for (int i=0; i<derivs.size(); ++i)
        derivs(i) = constraints[i]->transformDerivative(transformedParamters(i), derivs(i));

      return derivs;
    }

    calibration_adapter_ptr_t CalibrationAdapter::optionPricerAdapter( const opt_ptr_t& optionPricer,
                                                                       dbl_t strike,
                                                                       dbl_t forward,
                                                                       dbl_t expiry )
    {
      return std::make_shared<impl::OptionPricerAdapter>( optionPricer, strike, forward, expiry );
    }
  }
}