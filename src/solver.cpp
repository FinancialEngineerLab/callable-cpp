#include "solver.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "util.hpp"
#include "interpolation_builder.hpp"
#include "payoff.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct OneDimParabolicValuationPDESolver : public OneDimParabolicPDESolver
      {
        OneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                          const beagle::real_2d_function_ptr_t& diffusion) :
          m_Convection(convection),
          m_Diffusion(diffusion)
        { }
        virtual ~OneDimParabolicValuationPDESolver(void)
        { }
      public:
        virtual void evolve(double end,
                            double timeStep,
                            const beagle::dbl_vec_t& stateVariables,
                            const beagle::dbl_vec_t& lbc,
                            const beagle::dbl_vec_t& ubc,
                            beagle::dbl_vec_t& initialCondition) const override
        {
          // Assuming the state variable grid is uniform
          int numUnderlyingSteps = stateVariables.size();
          double deltaX = (stateVariables.back() - stateVariables.front()) / (numUnderlyingSteps - 1);

          // Now, perform induction
          beagle::dbl_vec_t diag(numUnderlyingSteps);
          beagle::dbl_vec_t lower(numUnderlyingSteps);
          beagle::dbl_vec_t upper(numUnderlyingSteps);

          double dTdX = timeStep / deltaX;
          double dTdXdX = dTdX / deltaX;

          for (int j=0; j<numUnderlyingSteps; ++j)
          {
            double convection = m_Convection->value(end, stateVariables[j]);
            double diffusion = m_Diffusion->value(end, stateVariables[j]);

            diag[j]  = 1. - 2. * diffusion * dTdXdX;
            upper[j] =   .5 * convection * dTdX + diffusion * dTdXdX;
            lower[j] = - .5 * convection * dTdX + diffusion * dTdXdX;
          }

          initialCondition[0] -= lower[0] * lbc[0];
          //diag[0]             -= lower[0] * lbc[2] / lbc[1];
          //upper[0]            -= lower[0] * lbc[3] / lbc[1];

          initialCondition[numUnderlyingSteps-1] -= upper[numUnderlyingSteps-1] * ubc[0];
          //diag[numUnderlyingSteps-1]             -= upper[numUnderlyingSteps-1] * ubc[2] / ubc[1];
          //upper[numUnderlyingSteps-1]            -= upper[numUnderlyingSteps-1] * ubc[3] / ubc[1];

          // for (auto price : prices)
          //   out << price << " ";
          // out << std::endl;

          beagle::util::tridiagonalSolve( initialCondition, diag, upper, lower );
        }
      private:
        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
        beagle::real_function_ptr_t m_Rate;
      };

      struct OneDimFokkerPlanckPDESolver : public OneDimParabolicPDESolver
      {
        OneDimFokkerPlanckPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                    const beagle::real_2d_function_ptr_t& diffusion) :
          m_Convection(convection),
          m_Diffusion(diffusion)
        { }
        virtual ~OneDimFokkerPlanckPDESolver(void)
        { }
      public:
        virtual void evolve(double end,
                            double timeStep,
                            const beagle::dbl_vec_t& stateVariables,
                            const beagle::dbl_vec_t& lbc,
                            const beagle::dbl_vec_t& ubc,
                            beagle::dbl_vec_t& initialCondition) const override
        { }
      private:
        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
      };
    }

    OneDimParabolicPDESolver::OneDimParabolicPDESolver(void)
    { }

    OneDimParabolicPDESolver::~OneDimParabolicPDESolver(void)
    { }

    beagle::parabolic_pde_solver_ptr_t
    OneDimParabolicPDESolver::formOneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                                    const beagle::real_2d_function_ptr_t& diffusion)
    {
      return std::make_shared<impl::OneDimParabolicValuationPDESolver>(convection, diffusion);
    }

    beagle::parabolic_pde_solver_ptr_t
    OneDimParabolicPDESolver::formOneDimFokkerPlanckPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                              const beagle::real_2d_function_ptr_t& diffusion)
    {
      return std::make_shared<impl::OneDimFokkerPlanckPDESolver>(convection, diffusion);
    }
  }
}
