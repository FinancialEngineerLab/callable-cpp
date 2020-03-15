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
        virtual beagle::dbl_t dividendAmount( beagle::dbl_t spot,
                                       beagle::dbl_t dividend ) const override
        {
          return spot > dividend ? dividend : spot;
        }
      };

      struct SurvivorDividendPolicy : public DividendPolicy
      {
        SurvivorDividendPolicy( void )
        { }
        virtual ~SurvivorDividendPolicy( void )
        { }
      public:
        virtual beagle::dbl_t dividendAmount( beagle::dbl_t spot,
                                       beagle::dbl_t dividend ) const override
        {
          return spot > dividend ? dividend : 0.;
        }
      };
    }

    DividendPolicy::DividendPolicy( void )
    { }

    DividendPolicy::~DividendPolicy( void )
    { }


    beagle::dbl_t
    DividendPolicy::exDividendStockPrice( beagle::dbl_t spot,
                                          beagle::dbl_t dividend ) const
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
