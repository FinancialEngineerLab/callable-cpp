#pragma once

#include "fwd_decl.hpp"

namespace beagle
{
  namespace product
  {
    namespace option
    {
      struct Payoff
      {
        virtual ~Payoff( void ) = default;
      public:
        virtual double intrinsicValue( double spot,
                                       double strike ) const = 0;
        virtual bool isCall( void ) const = 0;
        virtual bool isPut( void ) const = 0;
        virtual bool isDigital( void ) const = 0;
      public:
        static beagle::payoff_ptr_t call( void );
        static beagle::payoff_ptr_t put( void );
        static beagle::payoff_ptr_t digitalCall( void );
        static beagle::payoff_ptr_t digitalPut( void );
      };
    }
  }
}
