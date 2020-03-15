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
          CallPayoff( void ) :
            Payoff()
          { }
          virtual ~CallPayoff( void )
          { }
        public:
          virtual beagle::dbl_t intrinsicValue( beagle::dbl_t spot,
                                         beagle::dbl_t strike ) const override
          {
            return std::max<beagle::dbl_t>( spot - strike, 0. );
          }
          virtual bool isCall( void ) const override
          {
            return true;
          }
          virtual bool isPut( void ) const override
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
          virtual beagle::dbl_t intrinsicValue( beagle::dbl_t spot,
                                         beagle::dbl_t strike ) const override
          {
            return std::max<beagle::dbl_t>( strike - spot, 0. );
          }
          virtual bool isCall( void ) const override
          {
            return false;
          }
          virtual bool isPut( void ) const override
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
}
