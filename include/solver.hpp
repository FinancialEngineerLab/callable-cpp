#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct OneDimParabolicPDESolver
    {
      virtual ~OneDimParabolicPDESolver(void) = default;
    public:
      virtual void evolve(double end,
                          double timeStep,
                          const beagle::dbl_vec_t& stateVariables,
                          const beagle::dbl_vec_t& lowerBoundaryCondition,
                          const beagle::dbl_vec_t& upperBoundaryCondition,
                          beagle::dbl_vec_t& initialCondition) const = 0;
    public:
      static beagle::parabolic_pde_solver_ptr_t formOneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& drift,
                                                                                      const beagle::real_2d_function_ptr_t& vol,
                                                                                      const beagle::real_2d_function_ptr_t& rate,
                                                                                      const beagle::real_2d_function_ptr_t& source);
      static beagle::parabolic_pde_solver_ptr_t formOneDimFokkerPlanckPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                                                const beagle::real_2d_function_ptr_t& diffusion,
                                                                                const beagle::real_2d_function_ptr_t& rate);
      static beagle::parabolic_pde_solver_ptr_t formDupirePDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                                    const beagle::real_2d_function_ptr_t& diffusion);
    };
  }

}



#endif // !SOLVER_HPP
