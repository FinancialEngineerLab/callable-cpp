#ifndef ONE_DIM_PDE_PRICER_HPP
#define ONE_DIM_PDE_PRICER_HPP

#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "util.hpp"
#include "dividend_policy.hpp"
#include "real_function.hpp"
#include "interpolation_builder.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct OneDimensionalPDEOptionPricer : public Pricer,
                                             public beagle::valuation::mixins::FiniteDifference
      {
        explicit OneDimensionalPDEOptionPricer( const FiniteDifferenceDetails& fdDetails) :
          m_FDDetails(fdDetails)
        { }
        virtual ~OneDimensionalPDEOptionPricer( void )
        { }
      public:
        virtual const FiniteDifferenceDetails& finiteDifferenceDetails(void) const
        {
          return m_FDDetails;
        }
      private:
        FiniteDifferenceDetails m_FDDetails;
      };
    }
  }
}

#endif