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
          FixedCouponBond(double expiry,
                          double coupon,
                          int frequency)
          {
            m_Cashflows.clear();

            double accrual = 1. / frequency;
            double couponAmount = standardFaceValue() * coupon * accrual;

            int numPayments = std::ceil(expiry * frequency);
            for (int i=numPayments; i>1; --i)
            {
              m_Cashflows.emplace_back(expiry - (i-1) * accrual, couponAmount);
            }

            m_Cashflows.emplace_back(expiry, standardFaceValue() + couponAmount);
          }
          virtual ~FixedCouponBond( void )
          { }
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("FixedCouponBond");
            return ss;
          }
          virtual const beagle::bond_cashflows_t& cashflows(void) const override
          {
            return m_Cashflows;
          }
        private:
          beagle::bond_cashflows_t m_Cashflows;
        };
      }

      Bond::Bond(void)
      { }

      Bond::~Bond( void )
      { }

      beagle::product_ptr_t
      Bond::createFixedCouponBond(double expiry,
                                  double coupon,
                                  int frequency)
      {
        return std::make_shared<impl::FixedCouponBond>(expiry, coupon, frequency);
      }
    }
  }
}
