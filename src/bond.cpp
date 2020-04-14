#include "bond.hpp"

namespace beagle
{
  namespace product
  {
    namespace bond
    {
      namespace impl
      {
        struct ZeroCouponBond : public Bond
        {
          ZeroCouponBond(double expiry)
          {
            m_Couponflows.clear();
            m_Notionalflows.clear();
            m_Notionalflows.emplace_back(expiry, standardFaceValue());
          }
          virtual ~ZeroCouponBond( void ) = default;
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("ZeroCouponBond");
            return ss;
          }
          virtual const beagle::coupon_flows_t& couponFlows(void) const override
          {
            return m_Couponflows;
          }
          virtual const beagle::notional_flows_t& notionalFlows(void) const override
          {
            return m_Notionalflows;
          }
        private:
          beagle::coupon_flows_t m_Couponflows;
          beagle::notional_flows_t m_Notionalflows;
        };

        struct FixedCouponBond : public Bond
        {
          FixedCouponBond(double expiry,
                          double coupon,
                          int frequency)
          {
            m_Couponflows.clear();
            m_Notionalflows.clear();

            double accrual = 1. / frequency;
            double couponAmount = standardFaceValue() * coupon * accrual;

            int numPayments = static_cast<int>(std::ceil(expiry * frequency));
            for (int i=0; i<numPayments; ++i)
            {
              m_Couponflows.emplace_back(expiry - (numPayments-i-1) * accrual, couponAmount);
            }

            m_Notionalflows.emplace_back(expiry, standardFaceValue());
          }
          virtual ~FixedCouponBond( void ) = default;
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("FixedCouponBond");
            return ss;
          }
          virtual const beagle::coupon_flows_t& couponFlows(void) const override
          {
            return m_Couponflows;
          }
          virtual const beagle::notional_flows_t& notionalFlows(void) const override
          {
            return m_Notionalflows;
          }
        private:
          beagle::coupon_flows_t m_Couponflows;
          beagle::notional_flows_t m_Notionalflows;
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
          virtual ~ConvertibleBond( void ) = default;
        public:
          virtual const std::string& name(void) const override
          {
            static std::string ss("ConvertibleBond");
            return ss;
          }
          virtual const beagle::coupon_flows_t& couponFlows(void) const override
          {
            auto pB = dynamic_cast<beagle::product::mixins::Bond*>(m_UnderlyingBond.get());
            return pB->couponFlows();
          }
          virtual const beagle::notional_flows_t& notionalFlows(void) const override
          {
            auto pB = dynamic_cast<beagle::product::mixins::Bond*>(m_UnderlyingBond.get());
            return pB->notionalFlows();
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
