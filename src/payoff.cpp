#include "payoff.hpp"

namespace beagle
{
  namespace product
  {
    namespace option
    {
      namespace impl
      {
        struct CallPayoff : public Payoff
        {
        public:
          virtual double intrinsicValue( double spot,
                                         double strike ) const override
          {
            return std::max( spot - strike, 0. );
          }
          virtual bool isCall( void ) const override
          {
            return true;
          }
          virtual bool isPut( void ) const override
          {
            return false;
          }
          virtual bool isDigital( void ) const override
          {
            return false;
          }
        };

        struct PutPayoff : public Payoff
        {
        public:
          virtual double intrinsicValue( double spot,
                                         double strike ) const override
          {
            return std::max( strike - spot, 0. );
          }
          virtual bool isCall( void ) const override
          {
            return false;
          }
          virtual bool isPut( void ) const override
          {
            return true;
          }
          virtual bool isDigital( void ) const override
          {
            return false;
          }
        };

        struct DigitalCallPayoff : public Payoff
        {
        public:
          virtual double intrinsicValue( double spot,
                                         double strike ) const override
          {
            return spot > strike ? 1. : 0.;
          }
          virtual bool isCall( void ) const override
          {
            return true;
          }
          virtual bool isPut( void ) const override
          {
            return false;
          }
          virtual bool isDigital( void ) const override
          {
            return true;
          }
        };

        struct DigitalPutPayoff : public Payoff
        {
        public:
          virtual double intrinsicValue( double spot,
                                         double strike ) const override
          {
            return spot < strike ? 1. : 0.;
          }
          virtual bool isCall( void ) const override
          {
            return false;
          }
          virtual bool isPut( void ) const override
          {
            return true;
          }
          virtual bool isDigital( void ) const override
          {
            return true;
          }
        };
      }

      beagle::payoff_ptr_t
      Payoff::call( void )
      {
        return std::make_shared<impl::CallPayoff>();
      }

      beagle::payoff_ptr_t
      Payoff::put( void )
      {
        return std::make_shared<impl::PutPayoff>();
      }

      beagle::payoff_ptr_t
      Payoff::digitalCall( void )
      {
        return std::make_shared<impl::DigitalCallPayoff>();
      }

      beagle::payoff_ptr_t
      Payoff::digitalPut( void )
      {
        return std::make_shared<impl::DigitalPutPayoff>();
      }
    }
  }
}
