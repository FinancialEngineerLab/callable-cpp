#ifndef OPTION_HPP
#define OPTION_HPP

#include <memory>
#include "payoff.hpp"

namespace beagle
{
  struct Option;
  using option_ptr_t = std::shared_ptr<Option>;

  struct Option
  {
    Option( double expiry,
            double strike,
            const payoff_ptr_t& payoff );
    virtual ~Option( void );
  public:
    virtual double strike( void ) const;
    virtual double expiry( void ) const;
    virtual const payoff_ptr_t& payoff( void ) const;
  private:
    double m_Expiry;
    double m_Strike;
    payoff_ptr_t m_Payoff;

  public:
    static option_ptr_t createEuropeanOption( double expiry,
                                              double strike,
                                              const payoff_ptr_t& payoff );
    static option_ptr_t createAmericanOption( double expiry,
                                              double strike,
                                              const payoff_ptr_t& payoff );
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

#endif