#include "dividend_policy.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct LiquidatorDividendPolicy : public DividendPolicy
      {
        LiquidatorDividendPolicy( void )
        { }
        virtual ~LiquidatorDividendPolicy( void )
        { }
      public:
        virtual double exDividendStockPrice( double spot,
                                             double dividend ) const override
        {
          return spot > dividend ? spot - dividend : 0.;
        }
      };

      struct SurvivorDividendPolicy : public DividendPolicy
      {
        SurvivorDividendPolicy( void )
        { }
        virtual ~SurvivorDividendPolicy( void )
        { }
      public:
        virtual double exDividendStockPrice( double spot,
                                             double dividend ) const override
        {
          return spot > dividend ? spot - dividend : spot;
        }
      };
    }

    DividendPolicy::DividendPolicy( void )
    { }

    DividendPolicy::~DividendPolicy( void )
    { }

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