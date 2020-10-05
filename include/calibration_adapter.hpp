#pragma once

#include "fwd_decl.hpp"
#include "calibration_constraint.hpp"
#include "pricer.hpp"

namespace beagle
{
  namespace calibration
  {
    struct CalibrationAdapter
    {
      virtual ~CalibrationAdapter( void ) = default;
    public:
      virtual dbl_vec_t values( const dbl_vec_t& parameters ) const = 0;
      virtual dbl_mat_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const = 0;

      virtual dbl_vec_t values( const dbl_vec_t& transformedParamters,
                                const calibration_bound_constraint_coll_t& constraints ) const;
      virtual dbl_mat_t derivativeWithRespectToTransformedParameters(
                     const dbl_vec_t& transformedParamters,
                     const calibration_bound_constraint_coll_t& constraints ) const;
    };

    namespace util
    {
      dbl_vec_t getOriginalParameters( const dbl_vec_t& transformedParamters,
                                       const calibration_bound_constraint_coll_t& constraints );
      dbl_vec_t getTransformedParameters( const dbl_vec_t& parameters,
                                          const calibration_bound_constraint_coll_t& constraints );
      dbl_vec_t liveValues( const calibration_adapter_ptr_t& adapter,
                            const dbl_vec_t& parameters );
      dbl_vec_t liveValues( const calibration_adapter_ptr_t& adapter,
                            const dbl_vec_t& transformedParamters,
                            const calibration_bound_constraint_coll_t& constraints );
      dbl_mat_t jacobianMatrix( const calibration_adapter_ptr_t& adapter,
                                const dbl_vec_t& parameters );
      dbl_mat_t jacobianMatrix( const calibration_adapter_ptr_t& adapter,
                                const dbl_vec_t& transformedParamters,
                                const calibration_bound_constraint_coll_t& constraints );
    }
  }
}

