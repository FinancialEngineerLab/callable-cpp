#ifndef CALIBRATION_ADAPTER_HPP
#define CALIBRATION_ADAPTER_HPP

#include "typedefs.hpp"
#include "calibration_constraint.hpp"
#include "option_pricer.hpp"

namespace sabr_test
{
  namespace calibration
  {
    struct CalibrationAdapter;
    using calibration_adapter_ptr_t = std::shared_ptr<CalibrationAdapter>;
    using calibration_adapter_coll_t = std::vector<calibration_adapter_ptr_t>;
    using opt_ptr_t = option::opt_ptr_t;

    struct CalibrationAdapter
    {
      virtual ~CalibrationAdapter( void );
    public:
      virtual dbl_t value( const dbl_vec_t& parameters ) const = 0;
      virtual dbl_vec_t derivativeWithRespectToParameters( const dbl_vec_t& parameters ) const = 0;

      virtual dbl_t value( const dbl_vec_t& transformedParamters,
                           const calibration_bound_constraint_coll_t& constraints ) const;
      virtual dbl_vec_t derivativeWithRespectToTransformedParameters( 
                     const dbl_vec_t& transformedParamters,
                     const calibration_bound_constraint_coll_t& constraints ) const;
    public:
      static calibration_adapter_ptr_t optionPricerAdapter( const opt_ptr_t& optionPricer,
                                                            dbl_t strike,
                                                            dbl_t forward,
                                                            dbl_t expiry );
    };

    namespace util
    {
      dbl_vec_t getOriginalParameters( const dbl_vec_t& transformedParamters,
                                       const calibration_bound_constraint_coll_t& constraints );
      dbl_vec_t getTransformedParameters( const dbl_vec_t& parameters,
                                          const calibration_bound_constraint_coll_t& constraints );
      dbl_vec_t liveValues( const calibration_adapter_coll_t& adapters,
                            const dbl_vec_t& parameters );
      dbl_vec_t liveValues( const calibration_adapter_coll_t& adapters,
                            const dbl_vec_t& transformedParamters,
                            const calibration_bound_constraint_coll_t& constraints );
      dbl_mat_t jacobianMatrix( const calibration_adapter_coll_t& adapters,
                                const dbl_vec_t& parameters );
      dbl_mat_t jacobianMatrix( const calibration_adapter_coll_t& adapters,
                                const dbl_vec_t& transformedParamters,
                                const calibration_bound_constraint_coll_t& constraints );
    }
  }
}


#endif