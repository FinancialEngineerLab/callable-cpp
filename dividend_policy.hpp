#ifndef DIVIDEND_POLICY_HPP
#define DIVIDEND_POLICY_HPP

#include <memory>

namespace beagle
{
  struct DividendPolicy;
  using dividend_policy_ptr_t = std::shared_ptr<DividendPolicy>;

  struct DividendPolicy
  {
    DividendPolicy( void );
    virtual ~DividendPolicy( void );
  public:
    virtual double exDividendStockPrice( double spot,
                                         double dividend ) const = 0;
  public:
    static dividend_policy_ptr_t liquidator( void );
    static dividend_policy_ptr_t survivor( void );
  };
}


#endif