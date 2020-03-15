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
        virtual beagle::dbl_t strike(void) const=0;
        virtual beagle::dbl_t expiry(void) const=0;
        virtual const beagle::payoff_ptr_t& payoff(void) const=0;
      };
    }
  }
}

#endif