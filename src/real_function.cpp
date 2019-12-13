#include "real_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace mixins
    {
      InterpolationParameters::~InterpolationParameters( void )
      { }
    }

    namespace impl
    {
      InterpolatedFunction::InterpolatedFunction( const beagle::dbl_vec_t& xValues,
                                                  const beagle::dbl_vec_t& yValues ) :
        m_XValues( xValues ),
        m_YValues( yValues )
      { }

      InterpolatedFunction::~InterpolatedFunction( void )
      { }

      const beagle::dbl_vec_t& InterpolatedFunction::xParameters( void ) const
      {
        return m_XValues;
      }

      const beagle::dbl_vec_t& InterpolatedFunction::yParameters( void ) const
      {
        return m_YValues;
      }

      struct ConstantFunction : public RealFunction
      {
        explicit ConstantFunction( double constant ) :
          m_Const( constant )
        { }
        virtual ~ConstantFunction( void )
        { }
      public:
        virtual double value( double arg ) const override
        {
          return m_Const;
        }
      private:
        double m_Const;
      };

      struct UnaryFunction : public RealFunction
      {
        UnaryFunction( const beagle::real_func_t& func ) :
          m_Func( func )
        { }
        virtual ~UnaryFunction( void )
        { }
      public:
        virtual double value( double arg ) const override
        {
          return m_Func(arg);
        }
      private:
        beagle::real_func_t m_Func;
      };

      struct LinearWithFlatExtrapolationInterpolatedFunction : public InterpolatedFunction
      {
        LinearWithFlatExtrapolationInterpolatedFunction( const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues ) :
          InterpolatedFunction( xValues, yValues )
        { }
        virtual ~LinearWithFlatExtrapolationInterpolatedFunction( void )
        { }
      public:
        virtual double value( double arg ) const override
        {
          const beagle::dbl_vec_t& xValues(xParameters());
          const beagle::dbl_vec_t& yValues(yParameters());

          if (arg < xValues.front())
            return yValues.front();
          else if (arg >= xValues.back())
            return yValues.back();
          else
          {
            auto it = xValues.cbegin() + 1;
            auto itEnd = xValues.cend();
            auto jt = yValues.cbegin() + 1;
            while (it != itEnd)
            {
              if (arg < *it)
                break;
              else
              {
                ++it;
                ++jt;
              }
            }

            double xLeft = *(it - 1);
            double xRight = *it;
            double yLeft = *(jt - 1);
            double yRight = *jt;
            return yLeft + (arg - xLeft) / (xRight - xLeft) * (yRight - yLeft);
          }
        }
      };

      struct NaturalCubicSplineWithFlatExtrapolationInterpolatedFunction : public InterpolatedFunction
      {
        NaturalCubicSplineWithFlatExtrapolationInterpolatedFunction( const beagle::dbl_vec_t& xValues,
                                                                     const beagle::dbl_vec_t& yValues ) :
          InterpolatedFunction(xValues, yValues)
        { }
        virtual ~NaturalCubicSplineWithFlatExtrapolationInterpolatedFunction( void )
        { }
      public:
        virtual double value( double arg ) const override
        {
          const beagle::dbl_vec_t& xValues(xParameters());
          const beagle::dbl_vec_t& yValues(yParameters());

          return 0.0;
        }
      };

      struct PiecewiseConstantRightInterpolatedFunction : public InterpolatedFunction
      {
        PiecewiseConstantRightInterpolatedFunction( const beagle::dbl_vec_t& xValues,
                                                    const beagle::dbl_vec_t& yValues ) :
          InterpolatedFunction(xValues, yValues)
        { }
        virtual ~PiecewiseConstantRightInterpolatedFunction( void )
        { }
      public:
        virtual double value( double arg ) const override
        {
          const beagle::dbl_vec_t& xValues(xParameters());
          const beagle::dbl_vec_t& yValues(yParameters());

          if (arg >= xValues.back())
            return yValues.back();
          else
          {
            auto it = xValues.cbegin();
            auto itEnd = xValues.cend();
            auto jt = yValues.cbegin();
            while (it != itEnd)
            {
              if (arg < *it)
                break;
              else
              {
                ++it;
                ++jt;
              }
            }

            return *jt;
          }
        }
      };
    }


    RealFunction::RealFunction( void )
    { }

    RealFunction::~RealFunction( void )
    { }

    beagle::real_function_ptr_t
    RealFunction::createConstantFunction( double constant )
    {
      return std::make_shared<impl::ConstantFunction>( constant );
    }

    beagle::real_function_ptr_t
    RealFunction::createUnaryFunction( const beagle::real_func_t& func )
    {
      return std::make_shared<impl::UnaryFunction>( func );
    }

    beagle::real_function_ptr_t
    RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues )
    {
      if (xValues.size() != yValues.size())
        throw(std::string("Mismatch in the size of interpolation parameters"));

      if (xValues.size() == 1U)
        return createConstantFunction( yValues[0] );
      else
        return std::make_shared<impl::LinearWithFlatExtrapolationInterpolatedFunction>( xValues,
                                                                                        yValues );
    }

    beagle::real_function_ptr_t
    RealFunction::createNaturalCubicSplineWithFlatExtrapolationInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues )
    {
      if (xValues.size() != yValues.size())
        throw(std::string("Mismatch in the size of interpolation parameters"));

      if (xValues.size() == 1U)
        return createConstantFunction( yValues[0] );
      else
        return std::make_shared<impl::NaturalCubicSplineWithFlatExtrapolationInterpolatedFunction>( xValues,
                                                                                                    yValues );
    }

    beagle::real_function_ptr_t
    RealFunction::createPiecewiseConstantRightInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues )
    {
      if (xValues.size() != yValues.size())
        throw(std::string("Mismatch in the size of interpolation parameters"));

      if (xValues.size() == 1U)
        return createConstantFunction( yValues[0] );
      else
        return std::make_shared<impl::PiecewiseConstantRightInterpolatedFunction>( xValues,
                                                                                        yValues );
    }
  }
}