#include "solver.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"
#include "util.hpp"
#include "interpolation_builder.hpp"

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
        virtual beagle::dbl_vec_t boundaryCoefficients( void ) const override
        {
          return beagle::dbl_vec_t{0., 1., 0., 0.};
        }
      };

      struct LinearBoundaryCondition : public BoundaryCondition
      {
        virtual ~LinearBoundaryCondition( void )
        { }
      public:
        virtual beagle::dbl_vec_t boundaryCoefficients( void ) const override
        {
          return beagle::dbl_vec_t{0., 1., -2., 1.};
        }
      };

      struct OneDimParabolicValuationPDESolver : public OneDimParabolicPDESolver
      {
        OneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                          const beagle::real_2d_function_ptr_t& diffusion,
                                          const beagle::real_function_ptr_t& rate) :
          m_Convection(convection),
          m_Diffusion(diffusion),
          m_Rate(rate)
        { }
        virtual ~OneDimParabolicValuationPDESolver(void)
        { }
      public:
        virtual void evolve(double start,
                            double end,
                            int numTimeSteps,
                            const beagle::dbl_vec_t& stateVariables,
                            const beagle::boundary_condition_ptr_t& lowerBoundaryCondition,
                            const beagle::boundary_condition_ptr_t& upperBoundaryCondition,
                            beagle::dbl_vec_t& initialCondition) const override
        {
          // Form time steps
          beagle::dbl_vec_t times;
          double expiry = end - start;
          times.resize(numTimeSteps + 1);
          for (int i=0; i<numTimeSteps+1; ++i)
            times[i] = start + i * expiry / numTimeSteps;

          // Assuming the state variable grid is uniform
          int numUnderlyingSteps = stateVariables.size();
          double deltaX = (stateVariables.back() - stateVariables.front()) / (numUnderlyingSteps - 1);

          // Now, perform induction
          beagle::dbl_vec_t diag(numUnderlyingSteps);
          beagle::dbl_vec_t lower(numUnderlyingSteps);
          beagle::dbl_vec_t upper(numUnderlyingSteps);

          for (int i=0; i<numTimeSteps-1; ++i)
          {
            double thisTime = times[i+1];
            double deltaT = thisTime - times[i];
            double dTdX = deltaT / deltaX;
            double dTdXdX = dTdX / deltaX;

            for (int j=0; j<numUnderlyingSteps; ++j)
            {
              double convection = m_Convection->value(thisTime, stateVariables[j]);
              double diffusion = m_Diffusion->value(thisTime, stateVariables[j]);
              double rate = m_Rate->value(thisTime);

              diag[j]  = 1. - rate * deltaT - 2. * diffusion * dTdXdX;
              upper[j] =   .5 * convection * dTdX + diffusion * dTdXdX;
              lower[j] = - .5 * convection * dTdX + diffusion * dTdXdX;
            }

            beagle::dbl_vec_t lbc = lowerBoundaryCondition->boundaryCoefficients();
            initialCondition[0] -= lower[0] * lbc[0] / lbc[1];
            diag[0]             -= lower[0] * lbc[2] / lbc[1];
            upper[0]            -= lower[0] * lbc[3] / lbc[1];

            beagle::dbl_vec_t ubc = upperBoundaryCondition->boundaryCoefficients();
            initialCondition[numUnderlyingSteps-1] -= upper[numUnderlyingSteps-1] * ubc[0] / ubc[1];
            diag[numUnderlyingSteps-1]             -= upper[numUnderlyingSteps-1] * ubc[2] / ubc[1];
            upper[numUnderlyingSteps-1]            -= upper[numUnderlyingSteps-1] * ubc[3] / ubc[1];

            // for (auto price : prices)
            //   out << price << " ";
            // out << std::endl;

            beagle::util::tridiagonalSolve( initialCondition, diag, upper, lower );
          }
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
        virtual void evolve(double start,
                            double end,
                            int numTimeSteps,
                            const beagle::dbl_vec_t& stateVariables,
                            const beagle::boundary_condition_ptr_t& lowerBoundaryCondition,
                            const beagle::boundary_condition_ptr_t& upperBoundaryCondition,
                            beagle::dbl_vec_t& initialCondition) const override
        { }
      private:
        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
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

    OneDimParabolicPDESolver::OneDimParabolicPDESolver(void)
    { }

    OneDimParabolicPDESolver::~OneDimParabolicPDESolver(void)
    { }

    beagle::parabolic_pde_solver_ptr_t
    OneDimParabolicPDESolver::formOneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                                    const beagle::real_2d_function_ptr_t& diffusion,
                                                                    const beagle::real_function_ptr_t& rate)
    {
      return std::make_shared<impl::OneDimParabolicValuationPDESolver>(convection, diffusion, rate);
    }

    beagle::parabolic_pde_solver_ptr_t
    OneDimParabolicPDESolver::formOneDimFokkerPlanckPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                              const beagle::real_2d_function_ptr_t& diffusion)
    {
      return std::make_shared<impl::OneDimFokkerPlanckPDESolver>(convection, diffusion);
    }
  }
}
