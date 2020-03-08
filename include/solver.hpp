#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "pricer.hpp"

namespace beagle
{
  namespace math
  {
    struct BoundaryCondition
    {
      BoundaryCondition( void );
      virtual ~BoundaryCondition( void );
    public:
      virtual beagle::dbl_vec_t boundaryCoefficients( void ) const = 0;
    public:
      static beagle::boundary_condition_ptr_t naturalBoundaryCondition( void );
      static beagle::boundary_condition_ptr_t linearBoundaryCondition( void );
    };

    struct OneDimParabolicPDESolver
    {
      OneDimParabolicPDESolver(double spot,
                               const beagle::real_2d_function_ptr_t& drift,
                               const beagle::real_2d_function_ptr_t& volatility,
                               const beagle::real_function_ptr_t& rate,
                               const beagle::discrete_dividend_schedule_t& dividends,
                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
    public:
      void evolve(double start,
                  double end,
                  const beagle::boundary_condition_ptr_t& lowerBoundaryCondition,
                  const beagle::boundary_condition_ptr_t& upperBoundaryCondition,
                  beagle::real_function_ptr_t& initialValues) const;
    private:
      double m_Spot;
      beagle::real_2d_function_ptr_t m_Drift;
      beagle::real_2d_function_ptr_t m_Vol;
      beagle::real_function_ptr_t m_Rate;
      beagle::discrete_dividend_schedule_t m_Dividends;
      beagle::valuation::OneDimFiniteDifferenceSettings m_Settings;
    };
  }

}



#endif // !SOLVER_HPP
