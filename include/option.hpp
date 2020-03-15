#ifndef OPTION_HPP
#define OPTION_HPP

#include "fwd_decl.hpp"
#include "product.hpp"

namespace beagle
{
  namespace product
  {
    namespace option
    {
      struct Option : public Product,
                      public beagle::product::mixins::Option
      {
        Option(beagle::dbl_t expiry,
               beagle::dbl_t strike,
               const beagle::payoff_ptr_t& payoff);
        virtual ~Option(void);
      public:
        virtual const std::string& name(void) const = 0;
      public:
        virtual beagle::dbl_t strike(void) const;
        virtual beagle::dbl_t expiry(void) const;
        virtual const beagle::payoff_ptr_t& payoff(void) const;
      private:
        beagle::dbl_t m_Expiry;
        beagle::dbl_t m_Strike;
        beagle::payoff_ptr_t m_Payoff;

      public:
        static beagle::product_ptr_t createEuropeanOption(beagle::dbl_t expiry,
                                                          beagle::dbl_t strike,
                                                          const beagle::payoff_ptr_t& payoff);
        static beagle::product_ptr_t createAmericanOption(beagle::dbl_t expiry,
                                                          beagle::dbl_t strike,
                                                          const beagle::payoff_ptr_t& payoff);
      };

      namespace mixins
      {
        struct European
        {
          virtual ~European(void);
        };

        struct American
        {
          virtual ~American(void);
        };
      }
    }
  }
}

#endif