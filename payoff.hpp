#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <memory>
#include <string>

namespace beagle
{
  struct Payoff;
  using payoff_ptr_t = std::shared_ptr<Payoff>;

  struct Payoff
  {
    Payoff( void );
    virtual ~Payoff( void );
  public:
    virtual double intrinsicValue( double spot,
                                   double strike ) const = 0;
    virtual bool isCall( void ) const = 0;
    virtual bool isPut( void ) const = 0;
  public:
    static payoff_ptr_t call( void );
    static payoff_ptr_t put( void );
  };
}

#endif