#pragma once

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
        virtual const beagle::coupon_flows_t& couponFlows(void) const=0;
        virtual const beagle::notional_flows_t& notionalFlows(void) const=0;
      };
    }
  }
}
