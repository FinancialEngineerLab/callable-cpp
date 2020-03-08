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
      virtual beagle::dbl_vec_t lowerBoundaryCoefficients( void ) const = 0;
      virtual beagle::dbl_vec_t upperBoundaryCoefficients( void ) const = 0;
    public:
      static beagle::boundary_condition_ptr_t naturalBoundaryCondition( void );
      static beagle::boundary_condition_ptr_t linearBoundaryCondition( void );
    };

    struct OneDimParabolicPDESolver
    {
      OneDimParabolicPDESolver(double spot,
                               const beagle::real_2d_function_ptr_t& convection,
                               const beagle::real_2d_function_ptr_t& diffusion,
                               const beagle::real_function_ptr_t& rate,
                               const beagle::discrete_dividend_schedule_t& dividends,
                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
    public:
      beagle::real_function_ptr_t evolve(double start,
                                         double end,
                                         const beagle::boundary_condition_ptr_t& boundaryCondition,
                                         const beagle::real_function_ptr_t& initialCondition) const;
    private:
      double m_Spot;
      beagle::real_2d_function_ptr_t m_Convection;
      beagle::real_2d_function_ptr_t m_Diffusion;
      beagle::real_function_ptr_t m_Rate;
      beagle::discrete_dividend_schedule_t m_Dividends;
      beagle::valuation::OneDimFiniteDifferenceSettings m_Settings;
    };
  }

}



#endif // !SOLVER_HPP
