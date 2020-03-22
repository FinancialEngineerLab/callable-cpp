#ifndef BOND_HPP
#define BOND_HPP

#include "fwd_decl.hpp"
#include "product.hpp"

namespace beagle
{
  namespace product
  {
    namespace bond
    {
      struct Bond : public Product,
                    public beagle::product::mixins::Bond
      {
        Bond(const beagle::bond_cashflows_t& cashflows);
        virtual ~Bond(void);
      public:
        virtual const std::string& name(void) const=0;
      public:
        virtual const beagle::bond_cashflows_t& cashflows(void) const override;
      public:
        static beagle::product_ptr_t createFixedCouponBond(const beagle::bond_cashflows_t& cashflows);
      private:
        beagle::bond_cashflows_t m_Cashflows;

      };
    }
  }
}

#endif