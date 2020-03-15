#ifndef DIVIDEND_POLICY_HPP
#define DIVIDEND_POLICY_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace valuation
  {
    struct DividendPolicy
    {
      DividendPolicy( void );
      virtual ~DividendPolicy( void );
    public:
      virtual beagle::dbl_t exDividendStockPrice( beagle::dbl_t spot,
                                           beagle::dbl_t dividend ) const;
      virtual beagle::dbl_t dividendAmount( beagle::dbl_t spot,
                                     beagle::dbl_t dividend ) const = 0;
    public:
      static beagle::dividend_policy_ptr_t liquidator( void );
      static beagle::dividend_policy_ptr_t survivor( void );
    };
  }
}


#endif