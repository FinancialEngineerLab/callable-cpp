#include "calibration_constraint.hpp"
#include "util.hpp"

namespace beagle
{
  namespace calibration
  {
    namespace
    {
      beagle::dbl_t negativeInfiniteBoundary( void )
      {
        return -1e16;
      }

      beagle::dbl_t positiveInfiniteBoundary( void )
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
        virtual ~NoBoundCalibrationConstraintImpl( void )
        { }
      public:
        virtual beagle::dbl_t lowerBound( void ) const override
        {
          return negativeInfiniteBoundary();
        }
        virtual beagle::dbl_t upperBound( void ) const override
        {
          return positiveInfiniteBoundary();
        }
        virtual beagle::dbl_t transform( const beagle::dbl_t original ) const override
        {
          return original;
        }
        virtual beagle::dbl_t inverseTransform( const beagle::dbl_t transformed ) const override
        {
          return transformed;
        }
        virtual beagle::dbl_t transformDerivative( const beagle::dbl_t transformed,
                                            const beagle::dbl_t dOriginal ) const override
        {
          return dOriginal;
        }
      };

      struct LowerBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        explicit LowerBoundCalibrationConstraintImpl( const beagle::dbl_t lowerBound ) :
          m_LowerBound( lowerBound )
        { }
        virtual ~LowerBoundCalibrationConstraintImpl( void )
        { }
      public:
        virtual beagle::dbl_t lowerBound( void ) const override
        {
          return m_LowerBound;
        }
        virtual beagle::dbl_t upperBound( void ) const override
        {
          return positiveInfiniteBoundary();
        }
        virtual beagle::dbl_t transform( const beagle::dbl_t original ) const override
        {
          return std::log( original-m_LowerBound );
        }
        virtual beagle::dbl_t inverseTransform( const beagle::dbl_t transformed ) const override
        {
          return std::exp(transformed) + m_LowerBound;
        }
        virtual beagle::dbl_t transformDerivative( const beagle::dbl_t transformed,
                                            const beagle::dbl_t dOriginal ) const override
        {
          return dOriginal * std::exp(transformed);
        }
      private:
        const beagle::dbl_t m_LowerBound;
      };

      struct UpperBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        explicit UpperBoundCalibrationConstraintImpl( const beagle::dbl_t upperBound ) :
          m_UpperBound( upperBound )
        { }
        virtual ~UpperBoundCalibrationConstraintImpl( void )
        { }
      public:
        virtual beagle::dbl_t lowerBound( void ) const override
        {
          return negativeInfiniteBoundary();
        }
        virtual beagle::dbl_t upperBound( void ) const override
        {
          return m_UpperBound;
        }
        virtual beagle::dbl_t transform( const beagle::dbl_t original ) const override
        {
          return std::log( m_UpperBound-original );
        }
        virtual beagle::dbl_t inverseTransform( const beagle::dbl_t transformed ) const override
        {
          return m_UpperBound - std::exp(transformed);
        }
        virtual beagle::dbl_t transformDerivative( const beagle::dbl_t transformed,
                                            const beagle::dbl_t dOriginal ) const override
        {
          return -dOriginal * std::exp(transformed);
        }
      private:
        const beagle::dbl_t m_UpperBound;
      };

      struct TwoSidedBoundCalibrationConstraintImpl : public CalibrationBoundConstraint
      {
        TwoSidedBoundCalibrationConstraintImpl( const beagle::dbl_t lowerBound,
                                                const beagle::dbl_t upperBound ) :
          m_LowerBound( lowerBound ),
          m_UpperBound( upperBound ),
          m_UMinusL( upperBound-lowerBound )
        { }
        virtual ~TwoSidedBoundCalibrationConstraintImpl( void )
        { }
      public:
        virtual beagle::dbl_t lowerBound( void ) const override
        {
          return m_LowerBound;
        }
        virtual beagle::dbl_t upperBound( void ) const override
        {
          return m_UpperBound;
        }
        virtual beagle::dbl_t transform( const beagle::dbl_t original ) const override
        {
          return std::tan( util::pi() * ( (original-m_LowerBound)/m_UMinusL - 0.5 ) );
        }
        virtual beagle::dbl_t inverseTransform( const beagle::dbl_t transformed ) const override
        {
          return m_LowerBound + m_UMinusL * ( .5 * util::pi() + std::atan(transformed) ) / util::pi();
        }
        virtual beagle::dbl_t transformDerivative( const beagle::dbl_t transformed,
                                            const beagle::dbl_t dOriginal ) const override
        {
          return dOriginal * m_UMinusL / util::pi() / (1.0+transformed*transformed);
        }
      private:
        const beagle::dbl_t m_LowerBound;
        const beagle::dbl_t m_UpperBound;
        const beagle::dbl_t m_UMinusL;
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
    CalibrationBoundConstraint::lowerBoundCalibrationConstraint( const beagle::dbl_t lowerBound )
    {
      return std::make_shared<impl::LowerBoundCalibrationConstraintImpl>( lowerBound );
    }

    beagle::calibration_bound_constraint_ptr_t
    CalibrationBoundConstraint::upperBoundCalibrationConstraint( const beagle::dbl_t upperBound )
    {
      return std::make_shared<impl::UpperBoundCalibrationConstraintImpl>( upperBound );
    }

    beagle::calibration_bound_constraint_ptr_t
    CalibrationBoundConstraint::twoSidedBoundCalibrationConstraint( const beagle::dbl_t lowerBound,
                                                                    const beagle::dbl_t upperBound )
    {
      return std::make_shared<impl::TwoSidedBoundCalibrationConstraintImpl>( lowerBound, upperBound );
    }
  }

}
