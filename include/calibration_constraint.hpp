#ifndef CALIBRATION_CONSTRAINT_HPP
#define CALIBRATION_CONSTRAINT_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace calibration
  {
    struct CalibrationBoundConstraint
    {
    public:
      virtual ~CalibrationBoundConstraint( void );
      virtual beagle::dbl_t lowerBound( void ) const =0;
      virtual beagle::dbl_t upperBound( void ) const =0;
      virtual beagle::dbl_t transform( const beagle::dbl_t original ) const =0;
      virtual beagle::dbl_t inverseTransform( const beagle::dbl_t transformed ) const =0;
      virtual beagle::dbl_t transformDerivative( const beagle::dbl_t transformed,
                                                 const beagle::dbl_t dOriginal ) const =0;

    public:
      static beagle::calibration_bound_constraint_ptr_t noBoundCalibrationConstraint( void );
      static beagle::calibration_bound_constraint_ptr_t lowerBoundCalibrationConstraint( beagle::dbl_t lowerBound );
      static beagle::calibration_bound_constraint_ptr_t upperBoundCalibrationConstraint( beagle::dbl_t upperBound );
      static beagle::calibration_bound_constraint_ptr_t twoSidedBoundCalibrationConstraint( beagle::dbl_t lowerBound,
                                                                                            beagle::dbl_t upperBound );
    };
  }
}


#endif