#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace option
  {
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
      static beagle::payoff_ptr_t call( void );
      static beagle::payoff_ptr_t put( void );
    };
  }
}

#endif