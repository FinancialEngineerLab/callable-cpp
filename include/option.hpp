#ifndef OPTION_HPP
#define OPTION_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace option
  {
    struct Option
    {
      Option( double expiry,
              double strike,
              const beagle::payoff_ptr_t& payoff );
      virtual ~Option( void );
    public:
      virtual double strike( void ) const;
      virtual double expiry( void ) const;
      virtual const beagle::payoff_ptr_t& payoff( void ) const;
    private:
      double m_Expiry;
      double m_Strike;
      beagle::payoff_ptr_t m_Payoff;

    public:
      static beagle::option_ptr_t createEuropeanOption( double expiry,
                                                        double strike,
                                                        const beagle::payoff_ptr_t& payoff );
      static beagle::option_ptr_t createAmericanOption( double expiry,
                                                        double strike,
                                                        const beagle::payoff_ptr_t& payoff );
    };


    namespace mixins
    {
      struct European
      {
        virtual ~European( void );
      };

      struct American
      {
        virtual ~American( void );
      };
    }
  }
}

#endif