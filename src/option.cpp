#include "option.hpp"

namespace beagle
{
  namespace product
  {
    namespace option
    {
      namespace impl
      {
        struct EuropeanOption : public Option,
                                public beagle::product::option::mixins::European
        {
          EuropeanOption( beagle::dbl_t expiry,
                          beagle::dbl_t strike,
                          const beagle::payoff_ptr_t& payoff ) :
            Option( expiry, strike, payoff )
          { }
          virtual ~EuropeanOption( void )
          { }
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("EuropeanOption");
            return ss;
          }
        };

        struct AmericanOption : public Option,
                                public beagle::product::option::mixins::American
        {
          AmericanOption( beagle::dbl_t expiry,
                          beagle::dbl_t strike,
                          const beagle::payoff_ptr_t& payoff ) :
            Option( expiry, strike, payoff )
          { }
          virtual ~AmericanOption( void )
          { }
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("AmericanOption");
            return ss;
          }
        };
      }

      Option::Option( beagle::dbl_t expiry,
                      beagle::dbl_t strike,
                      const beagle::payoff_ptr_t& payoff ) :
        m_Expiry( expiry ),
        m_Strike( strike ),
        m_Payoff( payoff )
      { }

      Option::~Option( void )
      { }

      beagle::dbl_t Option::expiry( void ) const
      {
        return m_Expiry;
      }

      beagle::dbl_t Option::strike( void ) const
      {
        return m_Strike;
      }

      const beagle::payoff_ptr_t& Option::payoff( void ) const
      {
        return m_Payoff;
      }

      beagle::product_ptr_t
      Option::createEuropeanOption( beagle::dbl_t expiry,
                                    beagle::dbl_t strike,
                                    const beagle::payoff_ptr_t& payoff )
      {
        return std::make_shared<impl::EuropeanOption>( expiry, strike, payoff );
      }

      beagle::product_ptr_t
      Option::createAmericanOption( beagle::dbl_t expiry,
                                    beagle::dbl_t strike,
                                    const beagle::payoff_ptr_t& payoff )
      {
        return std::make_shared<impl::AmericanOption>( expiry, strike, payoff );
      }

      namespace mixins
      {
        European::~European( void ) { }
        American::~American( void ) { }
      }
    }
  }
}
