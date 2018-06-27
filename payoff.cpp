#include "payoff.hpp"

namespace beagle
{
  namespace option
  {
    namespace impl
    {
      struct CallPayoff : public Payoff
      {
        CallPayoff( void ) :
          Payoff()
        { }
        virtual ~CallPayoff( void )
        { }
      public:
        virtual double intrinsicValue( double spot,
                                       double strike ) const override
        {
          return std::max( spot - strike, 0. );
        }
        virtual bool isCall( void ) const
        {
          return true;
        }
        virtual bool isPut( void ) const
        {
          return false;
        }
      };

      struct PutPayoff : public Payoff
      {
        PutPayoff( void ) :
          Payoff()
        { }
        virtual ~PutPayoff( void )
        { }
      public:
        virtual double intrinsicValue( double spot,
                                       double strike ) const override
        {
          return std::max( strike - spot, 0. );
        }
        virtual bool isCall( void ) const
        {
          return false;
        }
        virtual bool isPut( void ) const
        {
          return true;
        }
      };
    }

    Payoff::Payoff( void )
    { }

    Payoff::~Payoff( void )
    { }

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
  }
}