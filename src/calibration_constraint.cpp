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
      public:
        virtual double lowerBound( void ) const override
        {
          return negativeInfiniteBoundary();
        }
        virtual double upperBound( void ) const override
        {
          return positiveInfiniteBoundary();
        }
        virtual double transform( const double original ) const override
        {
          return original;
        }
        virtual double inverseTransform( const double transformed ) const override
        {
          return transformed;
        }
        virtual double transformDerivative( const double transformed,
                                            const double dOriginal ) const override
        {
          return dOriginal;
        }
      };

      struct LowerBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        explicit LowerBoundCalibrationConstraintImpl(double lowerBound) :
          m_LowerBound( lowerBound )
        { }
      public:
        virtual double lowerBound( void ) const override
        {
          return m_LowerBound;
        }
        virtual double upperBound( void ) const override
        {
          return positiveInfiniteBoundary();
        }
        virtual double transform( const double original ) const override
        {
          return std::log( original-m_LowerBound );
        }
        virtual double inverseTransform( const double transformed ) const override
        {
          return std::exp(transformed) + m_LowerBound;
        }
        virtual double transformDerivative( const double transformed,
                                            const double dOriginal ) const override
        {
          return dOriginal * std::exp(transformed);
        }
      private:
        const double m_LowerBound;
      };

      struct UpperBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        explicit UpperBoundCalibrationConstraintImpl(double upperBound) :
          m_UpperBound( upperBound )
        { }
      public:
        virtual double lowerBound( void ) const override
        {
          return negativeInfiniteBoundary();
        }
        virtual double upperBound( void ) const override
        {
          return m_UpperBound;
        }
        virtual double transform( const double original ) const override
        {
          return std::log( m_UpperBound-original );
        }
        virtual double inverseTransform( const double transformed ) const override
        {
          return m_UpperBound - std::exp(transformed);
        }
        virtual double transformDerivative( const double transformed,
                                            const double dOriginal ) const override
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
      public:
        virtual double lowerBound( void ) const override
        {
          return m_LowerBound;
        }
        virtual double upperBound( void ) const override
        {
          return m_UpperBound;
        }
        virtual double transform( const double original ) const override
        {
          return std::tan( util::pi() * ( (original-m_LowerBound)/m_UMinusL - 0.5 ) );
        }
        virtual double inverseTransform( const double transformed ) const override
        {
          return m_LowerBound + m_UMinusL * ( .5 * util::pi() + std::atan(transformed) ) / util::pi();
        }
        virtual double transformDerivative( const double transformed,
                                            const double dOriginal ) const override
        {
          return dOriginal * m_UMinusL / util::pi() / (1.0+transformed*transformed);
        }
      private:
        const double m_LowerBound;
        const double m_UpperBound;
        const double m_UMinusL;
      };
    }

    const beagle::calibration_bound_constraint_ptr_t&
    CalibrationBoundConstraint::noBoundCalibrationConstraint( void )
    {
      static beagle::calibration_bound_constraint_ptr_t instance = std::make_shared<impl::NoBoundCalibrationConstraintImpl>();
      return instance;
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
