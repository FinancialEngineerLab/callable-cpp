#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "fwd_decl.hpp"

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
      OneDimParabolicPDESolver(void);
      virtual ~OneDimParabolicPDESolver(void);
    public:
      virtual void evolve(double start,
                          double end,
                          int numTimeSteps,
                          const beagle::dbl_vec_t& stateVariables,
                          const beagle::boundary_condition_ptr_t& lowerBoundaryCondition,
                          const beagle::boundary_condition_ptr_t& upperBoundaryCondition,
                          beagle::dbl_vec_t& initialCondition) const = 0;
    public:
      static beagle::parabolic_pde_solver_ptr_t formOneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                                                      const beagle::real_2d_function_ptr_t& diffusion,
                                                                                      const beagle::real_function_ptr_t& rate);
      static beagle::parabolic_pde_solver_ptr_t formOneDimFokkerPlanckPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                                                const beagle::real_2d_function_ptr_t& diffusion);
    private:
      beagle::real_2d_function_ptr_t m_Convection;
      beagle::real_2d_function_ptr_t m_Diffusion;
      beagle::real_function_ptr_t m_Rate;
    };
  }

}



#endif // !SOLVER_HPP
