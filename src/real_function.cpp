#include "real_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace mixins
    {
      InterpolationParameters::~InterpolationParameters( void )
      { }

      DividendSchedule::~DividendSchedule( void )
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
        explicit ConstantFunction( beagle::dbl_t constant ) :
          m_Const( constant )
        { }
        virtual ~ConstantFunction( void )
        { }
      public:
        virtual beagle::dbl_t value( beagle::dbl_t arg ) const override
        {
          return m_Const;
        }
      private:
        beagle::dbl_t m_Const;
      };

      struct CompositeFunction : public RealFunction
      {
        CompositeFunction( const beagle::real_function_ptr_t& f,
                           const beagle::real_function_ptr_t& g ) :
          m_F( f ),
          m_G( g )
        { }
        virtual ~CompositeFunction( void )
        { }
      public:
        virtual beagle::dbl_t value( beagle::dbl_t arg ) const override
        {
          return m_F->value( m_G->value(arg) );
        }
      private:
        beagle::real_function_ptr_t m_F;
        beagle::real_function_ptr_t m_G;
      };

      struct UnaryFunction : public RealFunction
      {
        explicit UnaryFunction(const beagle::real_func_t& func) :
          m_Func(func)
        { }
        virtual ~UnaryFunction(void)
        { }
      public:
        virtual beagle::dbl_t value(beagle::dbl_t arg) const override
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
        virtual beagle::dbl_t value( beagle::dbl_t arg ) const override
        {
          const beagle::dbl_vec_t& xValues(xParameters());
          const beagle::dbl_vec_t& yValues(yParameters());

          auto it = std::lower_bound(xValues.cbegin(),
                                     xValues.cend(),
                                     arg);

          if (it == xValues.cbegin())
            return yValues.front();
          else if (it == xValues.cend())
            return yValues.back();
          else
          {
            beagle::dbl_vec_t::difference_type diff = it - xValues.cbegin();
            auto jt = yValues.cbegin();
            std::advance(jt, diff);

            beagle::dbl_t xLeft = *(it - 1);
            beagle::dbl_t xRight = *it;
            beagle::dbl_t yLeft = *(jt - 1);
            beagle::dbl_t yRight = *jt;
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
        virtual beagle::dbl_t value( beagle::dbl_t arg ) const override
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
        virtual beagle::dbl_t value( beagle::dbl_t arg ) const override
        {
          const beagle::dbl_vec_t& xValues(xParameters());
          const beagle::dbl_vec_t& yValues(yParameters());

          auto it = std::lower_bound(xValues.cbegin(),
                                     xValues.cend(),
                                     arg);

          if (it == xValues.cend())
            return yValues.back();
          else
          {
            beagle::dbl_vec_t::difference_type diff = it - xValues.cbegin();
            auto jt = yValues.cbegin();
            std::advance(jt, diff);
            return *jt;
          }
        }
      };

      struct ContinuousForwardAssetPriceFunction : public RealFunction
      {
        ContinuousForwardAssetPriceFunction( beagle::dbl_t spot,
                                             const beagle::real_function_ptr_t& funding ) :
          m_Spot(spot),
          m_Funding(funding)
        { }
        virtual ~ContinuousForwardAssetPriceFunction( void )
        { }
      public:
        virtual beagle::dbl_t value( beagle::dbl_t arg ) const override
        {
          return m_Spot / m_Funding->value(arg);
        }
      private:
        beagle::dbl_t m_Spot;
        beagle::real_function_ptr_t m_Funding;
      };

      struct GeneralForwardAssetPriceFunction : public RealFunction,
                                                public mixins::DividendSchedule
      {
        GeneralForwardAssetPriceFunction( beagle::dbl_t spot,
                                          const beagle::real_function_ptr_t& funding,
                                          const beagle::discrete_dividend_schedule_t& dividends ) :
          m_Spot(spot),
          m_Funding(funding),
          m_Dividends(dividends)
        { }
        virtual ~GeneralForwardAssetPriceFunction( void )
        { }
      public:
        virtual beagle::dbl_t value( beagle::dbl_t arg ) const override
        {
          beagle::dbl_t factor = m_Funding->value(arg);
          beagle::dbl_t forward = m_Spot / factor;

          auto it = std::lower_bound(m_Dividends.cbegin(),
                                     m_Dividends.cend(),
                                     arg,
                                     [](const beagle::discrete_dividend_schedule_t::value_type& pair,
                                        beagle::dbl_t value)
                                     { return pair.first < value; });
          for (auto jt = m_Dividends.cbegin(); jt != it; ++jt)
          {
            beagle::dbl_t exDivTime = jt->first;
            beagle::dbl_t dividend = jt->second;
            beagle::dbl_t exDivFactor = m_Funding->value(exDivTime);
            forward -= dividend * exDivFactor / factor;
          }

          return forward;
        }
        virtual const beagle::discrete_dividend_schedule_t& dividendSchedule( void ) const override
        {
          return m_Dividends;
        }
      private:
        beagle::dbl_t m_Spot;
        beagle::real_function_ptr_t m_Funding;
        beagle::discrete_dividend_schedule_t m_Dividends;
      };
    }


    RealFunction::RealFunction( void )
    { }

    RealFunction::~RealFunction( void )
    { }

    beagle::real_function_ptr_t
    RealFunction::createConstantFunction( beagle::dbl_t constant )
    {
      return std::make_shared<impl::ConstantFunction>( constant );
    }

    beagle::real_function_ptr_t
    RealFunction::createUnaryFunction( const beagle::real_func_t& func )
    {
      return std::make_shared<impl::UnaryFunction>( func );
    }

    beagle::real_function_ptr_t
    RealFunction::createCompositeFunction( const beagle::real_function_ptr_t& f,
                                           const beagle::real_function_ptr_t& g )
    {
      return std::make_shared<impl::CompositeFunction>( f, g );
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

    beagle::real_function_ptr_t
    RealFunction::createContinuousForwardAssetPriceFunction( beagle::dbl_t spot,
                                                             const beagle::real_function_ptr_t& funding )
    {
      return std::make_shared<impl::ContinuousForwardAssetPriceFunction>( spot, funding );
    }

    beagle::real_function_ptr_t
    RealFunction::createGeneralForwardAssetPriceFunction( beagle::dbl_t spot,
                                                          const beagle::real_function_ptr_t& funding,
                                                          const beagle::discrete_dividend_schedule_t& dividends )
    {
      return std::make_shared<impl::GeneralForwardAssetPriceFunction>( spot, funding, dividends );
    }
  }
}