#include "dividend_policy.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct LiquidatorDividendPolicy : public DividendPolicy
      {
      public:
        virtual double dividendAmount( double spot,
                                       double dividend ) const override
        {
          return spot > dividend ? dividend : spot;
        }
      };

      struct SurvivorDividendPolicy : public DividendPolicy
      {
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

    const beagle::dividend_policy_ptr_t&
    DividendPolicy::liquidator( void )
    {
      static beagle::dividend_policy_ptr_t instance = std::make_shared<impl::LiquidatorDividendPolicy>();
      return instance;
    }

    const beagle::dividend_policy_ptr_t&
    DividendPolicy::survivor( void )
    {
      static beagle::dividend_policy_ptr_t instance = std::make_shared<impl::SurvivorDividendPolicy>();
      return instance;
    }
  }
}
