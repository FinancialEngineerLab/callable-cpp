#include "option.hpp"

namespace beagle
{
  namespace impl
  {
    struct EuropeanOption : public Option,
                            public beagle::mixins::European
    {
      EuropeanOption( double expiry,
                      double strike,
                      const payoff_ptr_t& payoff ) :
        Option( expiry, strike, payoff )
      { }
      virtual ~EuropeanOption( void )
      { }
    };

    struct AmericanOption : public Option,
                            public beagle::mixins::American
    {
      AmericanOption( double expiry,
                      double strike,
                      const payoff_ptr_t& payoff ) :
        Option( expiry, strike, payoff )
      { }
      virtual ~AmericanOption( void )
      { }
    };
  }

  Option::Option( double expiry,
                  double strike,
                  const payoff_ptr_t& payoff ) :
    m_Expiry( expiry ),
    m_Strike( strike ),
    m_Payoff( payoff )
  { }

  Option::~Option( void )
  { }

  double Option::expiry( void ) const
  {
    return m_Expiry;
  }

  double Option::strike( void ) const
  {
    return m_Strike;
  }

  const payoff_ptr_t& Option::payoff( void ) const
  {
    return m_Payoff;
  }

  option_ptr_t
  Option::createEuropeanOption( double expiry,
                                double strike,
                                const payoff_ptr_t& payoff )
  {
    return std::make_shared<impl::EuropeanOption>( expiry, strike, payoff );
  }

  option_ptr_t
  Option::createAmericanOption( double expiry,
                                double strike,
                                const payoff_ptr_t& payoff )
  {
    return std::make_shared<impl::AmericanOption>( expiry, strike, payoff );
  }

  namespace mixins
  {
    European::~European( void ) { }
    American::~American( void ) { }
  }
}