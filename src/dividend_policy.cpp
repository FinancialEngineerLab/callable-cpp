#include "dividend_policy.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct LiquidatorDividendPolicy : public DividendPolicy
      {
        virtual ~LiquidatorDividendPolicy( void ) = default;
      public:
        virtual double dividendAmount( double spot,
                                       double dividend ) const override
        {
          return spot > dividend ? dividend : spot;
        }
      };

      struct SurvivorDividendPolicy : public DividendPolicy
      {
        virtual ~SurvivorDividendPolicy( void ) = default;
      public:
        virtual double dividendAmount( double spot,
                                       double dividend ) const override
        {
          return spot > dividend ? dividend : 0.;
        }
      };
    }

    double
    DividendPolicy::exDividendStockPrice( double spot,
                                          double dividend ) const
    {
      return spot - dividendAmount( spot, dividend );
    }

    beagle::dividend_policy_ptr_t
    DividendPolicy::liquidator( void )
    {
      return std::make_shared<impl::LiquidatorDividendPolicy>();
    }

    beagle::dividend_policy_ptr_t
    DividendPolicy::survivor( void )
    {
      return std::make_shared<impl::SurvivorDividendPolicy>();
    }
  }
}
