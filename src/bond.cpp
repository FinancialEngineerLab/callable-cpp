#include "bond.hpp"

namespace beagle
{
  namespace product
  {
    namespace bond
    {
      namespace impl
      {
        struct FixedCouponBond : public Bond
        {
          FixedCouponBond(const beagle::bond_cashflows_t& cashflows) :
            Bond(cashflows)
          { }
          virtual ~FixedCouponBond( void )
          { }
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("FixedCouponBond");
            return ss;
          }
        };
      }

      Bond::Bond(const beagle::bond_cashflows_t& cashflows) :
        m_Cashflows(cashflows)
      { }

      Bond::~Bond( void )
      { }

      const beagle::bond_cashflows_t& Bond::cashflows(void) const
      {
        return m_Cashflows;
      }

      beagle::product_ptr_t
      Bond::createFixedCouponBond(const beagle::bond_cashflows_t& cashflows)
      {
        return std::make_shared<impl::FixedCouponBond>(cashflows);
      }
    }
  }
}
