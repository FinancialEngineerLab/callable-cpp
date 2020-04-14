#ifndef DIVIDEND_POLICY_HPP
#define DIVIDEND_POLICY_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace valuation
  {
    struct DividendPolicy
    {
      virtual ~DividendPolicy( void ) = default;
    public:
      virtual double exDividendStockPrice( double spot,
                                           double dividend ) const;
      virtual double dividendAmount( double spot,
                                     double dividend ) const = 0;
    public:
      static beagle::dividend_policy_ptr_t liquidator( void );
      static beagle::dividend_policy_ptr_t survivor( void );
    };
  }
}


#endif