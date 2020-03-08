#include "solver.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct NaturalBoundaryCondition : public BoundaryCondition
      {
        virtual ~NaturalBoundaryCondition( void )
        { }
      public:
        virtual beagle::dbl_vec_t boundaryCoefficients( void ) const
        {
          return beagle::dbl_vec_t{};
        }
      };
      
      struct LinearBoundaryCondition : public BoundaryCondition
      {
        virtual ~LinearBoundaryCondition( void )
        { }
      public:
        virtual beagle::dbl_vec_t boundaryCoefficients( void ) const
        {
          return beagle::dbl_vec_t{};
        }
      };
    }

    BoundaryCondition::BoundaryCondition( void )
    { }

    BoundaryCondition::~BoundaryCondition( void )
    { }

    beagle::boundary_condition_ptr_t
    BoundaryCondition::naturalBoundaryCondition( void )
    {
      return std::make_shared<impl::NaturalBoundaryCondition>();
    }

    beagle::boundary_condition_ptr_t
    BoundaryCondition::linearBoundaryCondition( void )
    {
      return std::make_shared<impl::LinearBoundaryCondition>();
    }

    OneDimParabolicPDESolver::OneDimParabolicPDESolver(double spot,
                                                       const beagle::real_2d_function_ptr_t& drift,
                                                       const beagle::real_2d_function_ptr_t& volatility,
                                                       const beagle::real_function_ptr_t& rate,
                                                       const beagle::discrete_dividend_schedule_t& dividends,
                                                       const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
      m_Spot(spot),
      m_Drift(drift),
      m_Vol(volatility),
      m_Rate(rate),
      m_Dividends(dividends),
      m_Settings(settings)
    { }

    void OneDimParabolicPDESolver::evolve(double start,
                                          double end,
                                          const beagle::boundary_condition_ptr_t& lowerBoundaryCondition,
                                          const beagle::boundary_condition_ptr_t& upperBoundaryCondition,
                                          beagle::real_function_ptr_t& initialValues) const
    {
      
    }
  }
}