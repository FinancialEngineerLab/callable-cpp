#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace product
  {
    namespace option
    {
      struct Payoff
      {
        Payoff( void );
        virtual ~Payoff( void );
      public:
        virtual beagle::dbl_t intrinsicValue( beagle::dbl_t spot,
                                       beagle::dbl_t strike ) const = 0;
        virtual bool isCall( void ) const = 0;
        virtual bool isPut( void ) const = 0;
      public:
        static beagle::payoff_ptr_t call( void );
        static beagle::payoff_ptr_t put( void );
      };
    }
  }
}

#endif