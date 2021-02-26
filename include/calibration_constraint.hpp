#pragma once

#include "fwd_decl.hpp"

namespace beagle
{
  namespace calibration
  {
    struct CalibrationBoundConstraint
    {
      virtual ~CalibrationBoundConstraint( void ) = default;
    public:
      virtual double lowerBound( void ) const =0;
      virtual double upperBound( void ) const =0;
      virtual double transform( const double original ) const =0;
      virtual double inverseTransform( const double transformed ) const =0;
      virtual double transformDerivative( const double transformed,
                                                 const double dOriginal ) const =0;

    public:
      static beagle::calibration_bound_constraint_ptr_t noBoundCalibrationConstraint( void );
      static beagle::calibration_bound_constraint_ptr_t lowerBoundCalibrationConstraint( double lowerBound );
      static beagle::calibration_bound_constraint_ptr_t upperBoundCalibrationConstraint( double upperBound );
      static beagle::calibration_bound_constraint_ptr_t twoSidedBoundCalibrationConstraint( double lowerBound,
                                                                                            double upperBound );
    };
  }
}
