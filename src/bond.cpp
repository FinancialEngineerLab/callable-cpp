#include "bond.hpp"

namespace beagle
{
  namespace product
  {
    namespace bond
    {
      namespace mixins
      {
        Callable::~Callable( void )
        { }

        Puttable::~Puttable( void )
        { }

        Convertible::~Convertible( void )
        { }
      }

      namespace impl
      {
        struct ZeroCouponBond : public Bond
        {
          ZeroCouponBond(double expiry)
          {
            m_Cashflows.clear();
            m_Cashflows.emplace_back(expiry, standardFaceValue());
          }
          virtual ~ZeroCouponBond( void )
          { }
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("ZeroCouponBond");
            return ss;
          }
          virtual const beagle::bond_cashflows_t& cashflows(void) const override
          {
            return m_Cashflows;
          }
        private:
          beagle::bond_cashflows_t m_Cashflows;
        };

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

        struct ConvertibleBond : public Bond,
                                 public beagle::product::bond::mixins::Callable,
                                 public beagle::product::bond::mixins::Puttable,
                                 public beagle::product::bond::mixins::Convertible
        {
          ConvertibleBond(const beagle::product_ptr_t& underlyingBond,
                          const beagle::real_function_ptr_t& conversionRatio,
                          const beagle::callable_schedule_t& callSchedule,
                          const beagle::puttable_schedule_t& putSchedule) :
            m_UnderlyingBond(underlyingBond),
            m_ConversonRatio(conversionRatio),
            m_CallSchedule(callSchedule),
            m_PutSchedule(putSchedule)
          { }
          virtual ~ConvertibleBond( void )
          { }
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("ConvertibleBond");
            return ss;
          }
          virtual const beagle::bond_cashflows_t& cashflows(void) const override
          {
            auto pB = dynamic_cast<beagle::product::mixins::Bond*>(m_UnderlyingBond.get());
            return pB->cashflows();
          }
          virtual const beagle::callable_schedule_t& callSchedule( void ) const override
          {
            return m_CallSchedule;
          }
          virtual const beagle::puttable_schedule_t& putSchedule( void ) const override
          {
            return m_PutSchedule;
          }
          virtual const beagle::real_function_ptr_t& conversionRatio( void ) const override
          {
            return m_ConversonRatio;
          }
        private:
          beagle::product_ptr_t m_UnderlyingBond;
          beagle::real_function_ptr_t m_ConversonRatio;
          beagle::callable_schedule_t m_CallSchedule;
          beagle::puttable_schedule_t m_PutSchedule;
        };
      }

      Bond::Bond(void)
      { }

      Bond::~Bond( void )
      { }

      beagle::product_ptr_t
      Bond::createZeroCouponBond(double expiry)
      {
        return std::make_shared<impl::ZeroCouponBond>(expiry);
      }

      beagle::product_ptr_t
      Bond::createFixedCouponBond(double expiry,
                                  double coupon,
                                  int frequency)
      {
        return std::make_shared<impl::FixedCouponBond>(expiry, coupon, frequency);
      }

      beagle::product_ptr_t
      Bond::createConvertibleBond(const beagle::product_ptr_t& underlyingBond,
                                  const beagle::real_function_ptr_t& conversionRatio,
                                  const beagle::callable_schedule_t& callSchedule,
                                  const beagle::puttable_schedule_t& putSchedule)
      {
        return std::make_shared<impl::ConvertibleBond>(underlyingBond, conversionRatio, callSchedule, putSchedule);
      }
    }
  }
}
