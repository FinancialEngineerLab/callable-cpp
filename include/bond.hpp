#pragma once

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
        virtual ~Bond(void) = default;
      public:
        virtual const std::string& name(void) const=0;
      public:
        virtual const beagle::coupon_flows_t& couponFlows(void) const=0;
        virtual const beagle::notional_flows_t& notionalFlows(void) const=0;
      public:
        static beagle::product_ptr_t createZeroCouponBond(double expiry);
        static beagle::product_ptr_t createFixedCouponBond(double expiry,
                                                           double coupon,
                                                           int frequency);
        static beagle::product_ptr_t createConvertibleBond(const beagle::product_ptr_t& underlyingBond,
                                                           const beagle::real_function_ptr_t& conversionRatio,
                                                           const beagle::callable_schedule_t& callSchedule,
                                                           const beagle::puttable_schedule_t& putSchedule);
      private:
        beagle::coupon_flows_t m_Cashflows;

      };

      namespace mixins
      {
        struct Callable
        {
          virtual ~Callable( void ) = default;
        public:
          virtual const beagle::callable_schedule_t& callSchedule( void ) const=0;
        };

        struct Puttable
        {
          virtual ~Puttable( void ) = default;
        public:
          virtual const beagle::puttable_schedule_t& putSchedule( void ) const=0;
        };

        struct Convertible
        {
          virtual ~Convertible( void ) = default;
        public:
          virtual const beagle::real_function_ptr_t& conversionRatio( void ) const=0;
        };
      }
    }
  }
}