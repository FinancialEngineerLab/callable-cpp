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
  }
}
