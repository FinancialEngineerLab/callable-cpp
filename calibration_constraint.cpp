#include "calibration_constraint.hpp"
#include "util.hpp"

namespace beagle
{
  namespace calibration
  {
    namespace
    {
      double negativeInfiniteBoundary( void )
      {
        return -1e16;
      }

      double positiveInfiniteBoundary( void )
      {
        return 1e16;
      }
    }

    namespace impl
    {
      struct NoBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        NoBoundCalibrationConstraintImpl( void )
        { }
        ~NoBoundCalibrationConstraintImpl( void )
        { }
      public:
        double lowerBound( void ) const
        {
          return negativeInfiniteBoundary();
        }
        double upperBound( void ) const
        {
          return positiveInfiniteBoundary();
        }
        double transform( const double original ) const
        {
          return original;
        }
        double inverseTransform( const double transformed ) const
        {
          return transformed;
        }
        double transformDerivative( const double transformed,
                                           const double dOriginal ) const
        {
          return dOriginal;
        }
      };

      struct LowerBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        LowerBoundCalibrationConstraintImpl( const double lowerBound ) :
          m_LowerBound( lowerBound )
        { }
        ~LowerBoundCalibrationConstraintImpl( void )
        { }
      public:
        double lowerBound( void ) const
        {
          return m_LowerBound;
        }
        double upperBound( void ) const
        {
          return positiveInfiniteBoundary();
        }
        double transform( const double original ) const
        {
          return std::log( original-m_LowerBound );
        }
        double inverseTransform( const double transformed ) const
        {
          return std::exp(transformed) + m_LowerBound;
        }
        double transformDerivative( const double transformed,
                                           const double dOriginal ) const
        {
          return dOriginal * std::exp(transformed);
        }
      private:
        const double m_LowerBound;
      };

      struct UpperBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        UpperBoundCalibrationConstraintImpl( const double upperBound ) :
          m_UpperBound( upperBound )
        { }
        ~UpperBoundCalibrationConstraintImpl( void )
        { }
      public:
        double lowerBound( void ) const
        {
          return negativeInfiniteBoundary();
        }
        double upperBound( void ) const
        {
          return m_UpperBound;
        }
        double transform( const double original ) const
        {
          return std::log( m_UpperBound-original );
        }
        double inverseTransform( const double transformed ) const
        {
          return m_UpperBound - std::exp(transformed);
        }
        double transformDerivative( const double transformed,
                                           const double dOriginal ) const
        {
          return -dOriginal * std::exp(transformed);
        }
      private:
        const double m_UpperBound;
      };

      struct TwoSidedBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        TwoSidedBoundCalibrationConstraintImpl( const double lowerBound,
                                                const double upperBound ) :
          m_LowerBound( lowerBound ),
          m_UpperBound( upperBound ),
          m_UMinusL( upperBound-lowerBound )
        { }
        ~TwoSidedBoundCalibrationConstraintImpl( void )
        { }
      public:
        double lowerBound( void ) const
        {
          return m_LowerBound;
        }
        double upperBound( void ) const
        {
          return m_UpperBound;
        }
        double transform( const double original ) const
        {
          return std::tan( util::pi() * ( (original-m_LowerBound)/m_UMinusL - 0.5 ) );
        }
        double inverseTransform( const double transformed ) const
        {
          return m_LowerBound + m_UMinusL * ( .5 * util::pi() + std::atan(transformed) ) / util::pi();
        }
        double transformDerivative( const double transformed,
                                           const double dOriginal ) const
        {
          return dOriginal * m_UMinusL / util::pi() / (1.0+transformed*transformed);
        }
      private:
        const double m_LowerBound;
        const double m_UpperBound;
        const double m_UMinusL;
      };
    }

    CalibrationBoundConstraint::~CalibrationBoundConstraint( void )
    { }

    beagle::calibration_bound_constraint_ptr_t
    CalibrationBoundConstraint::noBoundCalibrationConstraint( void )
    {
      return std::make_shared<impl::NoBoundCalibrationConstraintImpl>();
    }

    beagle::calibration_bound_constraint_ptr_t
    CalibrationBoundConstraint::lowerBoundCalibrationConstraint( const double lowerBound )
    {
      return std::make_shared<impl::LowerBoundCalibrationConstraintImpl>( lowerBound );
    }

    beagle::calibration_bound_constraint_ptr_t
    CalibrationBoundConstraint::upperBoundCalibrationConstraint( const double upperBound )
    {
      return std::make_shared<impl::UpperBoundCalibrationConstraintImpl>( upperBound );
    }

    beagle::calibration_bound_constraint_ptr_t
    CalibrationBoundConstraint::twoSidedBoundCalibrationConstraint( const double lowerBound,
                                                                    const double upperBound )
    {
      return std::make_shared<impl::TwoSidedBoundCalibrationConstraintImpl>( lowerBound, upperBound );
    }
  }

}