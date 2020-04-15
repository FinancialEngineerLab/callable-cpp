#include "real_function.hpp"
#include "dividend_policy.hpp"
#include "util.hpp"

#include <iostream>

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      InterpolatedFunction::InterpolatedFunction( const beagle::dbl_vec_t& xValues,
                                                  const beagle::dbl_vec_t& yValues ) :
        m_XValues( xValues ),
        m_YValues( yValues )
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
        virtual ~ConstantFunction( void ) = default;
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
        virtual ~CompositeFunction( void ) = default;
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
        virtual ~UnaryFunction(void) = default;
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
        virtual ~LinearWithFlatExtrapolationInterpolatedFunction( void ) = default;
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

      struct NaturalCubicSplineInterpolatedFunction : public InterpolatedFunction
      {
        NaturalCubicSplineInterpolatedFunction( const beagle::dbl_vec_t& xValues,
                                                const beagle::dbl_vec_t& yValues ) :
          InterpolatedFunction(xValues, yValues)
        {
          // By construction, the number of interpolation parameters is at least 3
          int sz = xValues.size();
          beagle::dbl_vec_t diag(sz-2, 0.);
          beagle::dbl_vec_t upper(sz-2, 0.);
          beagle::dbl_vec_t lower(sz-2, 0.);
          beagle::dbl_vec_t rhs(sz-2, 0.);
          
          for (int i=0; i<sz-2; ++i)
          {
            double h1 = xValues[i+1] - xValues[i];
            double h2 = xValues[i+2] - xValues[i+1];
            diag[i] = (h1 + h2) / 3.;
            upper[i] = h2 / 6.;
            lower[i] = h1 / 6.;
            rhs[i] = (yValues[i+2] - yValues[i+1]) / h2
                   - (yValues[i+1] - yValues[i]) / h1;
          }

          beagle::util::tridiagonalSolve( rhs, diag, upper, lower );

          m_SecondDerivs.resize(xValues.size(), 0.);
          std::copy(rhs.cbegin(),
                    rhs.cend(),
                    m_SecondDerivs.begin() + 1U);
        }
        virtual ~NaturalCubicSplineInterpolatedFunction( void ) = default;
      public:
        virtual double value( double arg ) const override
        {
          const beagle::dbl_vec_t& xValues(xParameters());
          const beagle::dbl_vec_t& yValues(yParameters());

          auto it = std::lower_bound(xValues.cbegin(),
                                     xValues.cend(),
                                     arg);

          double x1;
          double x2;
          double y1;
          double y2;
          double d1;
          double d2;
          if (it == xValues.cend())
          {
            x1 = *(it - 2U);
            x2 = *(it - 1U);

            auto jt = yValues.cend();
            y1 = *(jt - 2U);
            y2 = *(jt - 1U);

            auto kt = m_SecondDerivs.cend();
            d1 = *(kt - 2U);
            d2 = *(kt - 1U);
          }
          else if (it == xValues.cbegin())
          {
            x1 = *it;
            x2 = *(it + 1U);

            auto jt = yValues.cbegin();
            y1 = *jt;
            y2 = *(jt + 1U);

            auto kt = m_SecondDerivs.cbegin();
            d1 = *kt;
            d2 = *(kt + 1U);
          }
          else
          {
            x1 = *(it - 1U);
            x2 = *it;

            auto diff = std::distance(xValues.cbegin(), it);
            auto jt = yValues.cbegin();
            std::advance(jt, diff);
            y1 = *(jt - 1U);
            y2 = *jt;

            auto kt = m_SecondDerivs.cbegin();
            std::advance(kt, diff);
            d1 = *(kt - 1U);
            d2 = *kt;
          }

          double delta1 = x2 - arg;
          double delta2 = arg - x1;
          double h = x2 - x1;

          return std::pow(delta1, 3.) * d1 / 6. / h
               + std::pow(delta2, 3.) * d2 / 6. / h
               + delta1 * (y1 / h - h * d1 / 6.)
               + delta2 * (y2 / h - h * d2 / 6.);
        }
      private:
        beagle::dbl_vec_t m_SecondDerivs; 
      };

      struct PiecewiseConstantRightInterpolatedFunction : public InterpolatedFunction
      {
        PiecewiseConstantRightInterpolatedFunction( const beagle::dbl_vec_t& xValues,
                                                    const beagle::dbl_vec_t& yValues ) :
          InterpolatedFunction(xValues, yValues)
        { }
        virtual ~PiecewiseConstantRightInterpolatedFunction( void ) = default;
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
        virtual ~ContinuousForwardAssetPriceFunction( void ) = default;
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
        virtual ~GeneralForwardAssetPriceFunction( void ) = default;
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
    RealFunction::createNaturalCubicSplineInterpolatedFunction(
                                                         const beagle::dbl_vec_t& xValues,
                                                         const beagle::dbl_vec_t& yValues )
    {
      if (xValues.size() != yValues.size())
        throw(std::string("Mismatch in the size of interpolation parameters"));

      if (xValues.size() == 1U)
        return createConstantFunction( yValues[0] );
      else if (xValues.size() == 2U)
        return createLinearWithFlatExtrapolationInterpolatedFunction(xValues, yValues);
      else
        return std::make_shared<impl::NaturalCubicSplineInterpolatedFunction>( xValues,  yValues );
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