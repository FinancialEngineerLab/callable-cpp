#include "real_function.hpp"
#include "dividend_policy.hpp"

namespace beagle
{
  namespace math
  {
    namespace mixins
    {
      InterpolationParameters::~InterpolationParameters( void )
      { }

      ContainsDividends::~ContainsDividends( void )
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
        virtual double value( double arg ) const override
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
        virtual double value(double arg) const override
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
        ContinuousForwardAssetPriceFunction( double spot,
                                             const beagle::real_function_ptr_t& funding ) :
          m_Spot(spot),
          m_Funding(funding)
        { }
        virtual ~ContinuousForwardAssetPriceFunction( void )
        { }
      public:
        virtual double value( double arg ) const override
        {
          return m_Spot / m_Funding->value(arg);
        }
      private:
        double m_Spot;
        beagle::real_function_ptr_t m_Funding;
      };

      struct GeneralForwardAssetPriceFunction : public RealFunction,
                                                public mixins::ContainsDividends
      {
        GeneralForwardAssetPriceFunction( double spot,
                                          const beagle::real_function_ptr_t& funding,
                                          const beagle::dividend_schedule_t& dividends,
                                          const beagle::dividend_policy_ptr_t& policy ) :
          m_Spot(spot),
          m_Funding(funding),
          m_Dividends(dividends),
          m_Policy(policy)
        {
          m_Forwards.clear();

          double start = 1.;
          double forward = m_Spot;
          for (auto jt = m_Dividends.cbegin(); jt != m_Dividends.cend(); ++jt)
          {
            double end = m_Funding->value(std::get<0>(*jt));
            double rel = std::get<1>(*jt);
            double abs = std::get<2>(*jt);

            typename beagle::cum_ex_dividend_prices_t::value_type pair;

            forward *= start / end;
            pair.first = forward;

            forward *= (1. - rel);
            forward = m_Policy->exDividendStockPrice(forward, abs);
            pair.second = forward;

            start = end;
            m_Forwards.emplace_back(pair);
          }
        }
        virtual ~GeneralForwardAssetPriceFunction( void )
        { }
      public:
        virtual double value( double arg ) const override
        {
          auto it = std::lower_bound(m_Dividends.cbegin(),
                                     m_Dividends.cend(),
                                     arg,
                                     [](const beagle::dividend_schedule_t::value_type& dividend,
                                        double value)
                                     { return std::get<0>(dividend) < value; });

          if (it == m_Dividends.cbegin())
            return m_Spot / m_Funding->value(arg);
          else
          {
            --it;
            auto diff = std::distance(m_Dividends.cbegin(), it);
            auto jt = m_Forwards.cbegin();
            std::advance(jt, diff);
            return jt->second * m_Funding->value(std::get<0>(*it)) / m_Funding->value(arg);
          }
        }
        virtual const beagle::dividend_schedule_t& dividendSchedule( void ) const override
        {
          return m_Dividends;
        }
        virtual const beagle::dividend_policy_ptr_t& dividendPolicy( void ) const override
        {
          return m_Policy;
        }
        virtual const beagle::cum_ex_dividend_prices_t& cumAndExDividendForwards( void ) const override
        {
          return m_Forwards;
        }
      private:
        double m_Spot;
        beagle::real_function_ptr_t m_Funding;
        beagle::dividend_schedule_t m_Dividends;
        beagle::dividend_policy_ptr_t m_Policy;
        beagle::cum_ex_dividend_prices_t m_Forwards;
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
    RealFunction::createContinuousForwardAssetPriceFunction( double spot,
                                                             const beagle::real_function_ptr_t& funding )
    {
      return std::make_shared<impl::ContinuousForwardAssetPriceFunction>( spot, funding );
    }

    beagle::real_function_ptr_t
    RealFunction::createGeneralForwardAssetPriceFunction( double spot,
                                                          const beagle::real_function_ptr_t& funding,
                                                          const beagle::dividend_schedule_t& dividends,
                                                          const beagle::dividend_policy_ptr_t& policy )
    {
      if (dividends.empty())
        return RealFunction::createContinuousForwardAssetPriceFunction(spot, funding);

      return std::make_shared<impl::GeneralForwardAssetPriceFunction>( spot, funding, dividends, policy );
    }
  }
}