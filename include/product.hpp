#ifndef PRODUCT_HPP
#define PRODUCT_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace product
  {
    struct Product
    {
      virtual ~Product(void);
    public:
      virtual const std::string& name(void) const = 0;
    };

    namespace mixins
    {
      struct Option
      {
        virtual ~Option(void);
      public:
        virtual double strike(void) const=0;
        virtual double expiry(void) const=0;
        virtual const beagle::payoff_ptr_t& payoff(void) const=0;
      };

      struct Bond
      {
        virtual ~Bond(void);
      public:
        virtual double standardFaceValue(void) const;
        virtual const beagle::bond_cashflows_t& cashflows(void) const=0;
      };
    }
  }
}

#endif