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
        Option(double expiry,
               double strike,
               const beagle::payoff_ptr_t& payoff);
        virtual ~Option(void);
      public:
        virtual const std::string& name(void) const = 0;
      public:
        virtual double strike(void) const override;
        virtual double expiry(void) const override;
        virtual const beagle::payoff_ptr_t& payoff(void) const override;
      private:
        double m_Expiry;
        double m_Strike;
        beagle::payoff_ptr_t m_Payoff;
      public:
        static beagle::product_ptr_t createEuropeanOption(double expiry,
                                                          double strike,
                                                          const beagle::payoff_ptr_t& payoff);
        static beagle::product_ptr_t createAmericanOption(double expiry,
                                                          double strike,
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