#include "solver.hpp"
#include "pricer.hpp"
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
      struct OneDimParabolicPDEPricer : public beagle::valuation::Pricer
      {
        OneDimParabolicPDEPricer(const beagle::real_function_ptr_t& forward,
                                 const beagle::real_function_ptr_t& discounting,
                                 const beagle::real_2d_function_ptr_t& drift,
                                 const beagle::real_2d_function_ptr_t& volatility,
                                 const beagle::real_2d_function_ptr_t& rate,
                                 const beagle::real_2d_function_ptr_t& recovery,
                                 const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
          m_Forward(forward),
          m_Discounting(discounting),
          m_Vol(volatility),
          m_Settings(settings)
        {
          m_Convection = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              double localVol = volatility->value(time, spot);
                              double dft = drift->value(time, spot);
                              return dft - .5 * localVol * localVol;
                            } );
          m_Diffusion = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              double localVol = volatility->value(time, spot);
                              return .5 * localVol * localVol;
                            } );
          m_Rate = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              return rate->value(time, spot);
                            } );
          m_Source = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              return discounting->value(time) * recovery->value(time, spot);
                            } );}
        virtual ~OneDimParabolicPDEPricer( void ) = default;
      protected:
        const beagle::real_function_ptr_t& forwardCurve( void ) const
        {
          return m_Forward;
        }
        const beagle::real_function_ptr_t& discountCurve( void ) const
        {
          return m_Discounting;
        }
        const beagle::real_2d_function_ptr_t& volatilitySurface( void ) const
        {
          return m_Vol;
        }
        const beagle::valuation::OneDimFiniteDifferenceSettings& finiteDifferenceSettings(void) const
        {
          return m_Settings;
        }
        const beagle::real_2d_function_ptr_t& convectionCoefficient( void ) const
        {
          return m_Convection;
        }
        const beagle::real_2d_function_ptr_t& diffusionCoefficient( void ) const
        {
          return m_Diffusion;
        }
        const beagle::real_2d_function_ptr_t& rateCoefficient( void ) const
        {
          return m_Rate;
        }
        const beagle::real_2d_function_ptr_t& sourceTerm( void ) const
        {
          return m_Source;
        }
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        beagle::real_2d_function_ptr_t m_Vol;
        beagle::valuation::OneDimFiniteDifferenceSettings m_Settings;

        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
        beagle::real_2d_function_ptr_t m_Rate;
        beagle::real_2d_function_ptr_t m_Source;
      };

      struct OneDimParabolicValuationPDESolver : public OneDimParabolicPDESolver
      {
        OneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                          const beagle::real_2d_function_ptr_t& diffusion,
                                          const beagle::real_2d_function_ptr_t& rate,
                                          const beagle::real_2d_function_ptr_t& source) :
          m_Convection(convection),
          m_Diffusion(diffusion),
          m_Rate(rate),
          m_Source(source)
        { }
        virtual ~OneDimParabolicValuationPDESolver(void) = default;
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

          double dTdX = timeStep / deltaX;
          double dTdXdX = dTdX / deltaX;

          // Now, perform induction
          beagle::dbl_vec_t diag(numUnderlyingSteps);
          beagle::dbl_vec_t lower(numUnderlyingSteps);
          beagle::dbl_vec_t upper(numUnderlyingSteps);


          for (int j=0; j<numUnderlyingSteps; ++j)
          {
            double rate = m_Rate->value(end, stateVariables[j]);
            double convection = m_Convection->value(end, stateVariables[j]);
            double diffusion = m_Diffusion->value(end, stateVariables[j]);
            double source = m_Source->value(end, stateVariables[j]);

            diag[j]  = 1. - rate * timeStep - 2. * diffusion * dTdXdX;
            upper[j] =   .5 * convection * dTdX + diffusion * dTdXdX;
            lower[j] = - .5 * convection * dTdX + diffusion * dTdXdX;
            initialCondition[j] += source * timeStep;
          }

          initialCondition[0] -= lower[0] * lbc[0];
          //diag[0]             -= lower[0] * lbc[2] / lbc[1];
          //upper[0]            -= lower[0] * lbc[3] / lbc[1];

          initialCondition[numUnderlyingSteps-1] -= upper[numUnderlyingSteps-1] * ubc[0];
          //diag[numUnderlyingSteps-1]             -= upper[numUnderlyingSteps-1] * ubc[2] / ubc[1];
          //upper[numUnderlyingSteps-1]            -= upper[numUnderlyingSteps-1] * ubc[3] / ubc[1];

          beagle::util::tridiagonalSolve( initialCondition, diag, upper, lower );
        }
      private:
        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
        beagle::real_2d_function_ptr_t m_Rate;
        beagle::real_2d_function_ptr_t m_Source;
      };

      struct OneDimFokkerPlanckPDESolver : public OneDimParabolicPDESolver
      {
        OneDimFokkerPlanckPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                    const beagle::real_2d_function_ptr_t& diffusion,
                                    const beagle::real_2d_function_ptr_t& rate) :
          m_Convection(convection),
          m_Diffusion(diffusion),
          m_Rate(rate)
        { }
        virtual ~OneDimFokkerPlanckPDESolver(void) = default;
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
            double rate = .0; // m_Rate->value(end, stateVariables[j]);

            diag[j]  = 1.
                     + m_Rate->value(end, stateVariables[j]) * timeStep
                     + 2. * m_Diffusion->value(end, stateVariables[j]) * dTdXdX;
            if (j == 0)
              lower[j] = .0;
            else
              lower[j] = - .5 * m_Convection->value(end, stateVariables[j-1]) * dTdX
                         - m_Diffusion->value(end, stateVariables[j-1]) * dTdXdX;

            if (j == numUnderlyingSteps-1)
              upper[j] = .0;
            else
              upper[j] =   .5 * m_Convection->value(end, stateVariables[j+1]) * dTdX
                         - m_Diffusion->value(end, stateVariables[j+1]) * dTdXdX;
          }

          initialCondition[0] -= lower[0] * lbc[0];
          //diag[0]             -= lower[0] * lbc[2] / lbc[1];
          //upper[0]            -= lower[0] * lbc[3] / lbc[1];

          initialCondition[numUnderlyingSteps-1] -= upper[numUnderlyingSteps-1] * ubc[0];
          //diag[numUnderlyingSteps-1]             -= upper[numUnderlyingSteps-1] * ubc[2] / ubc[1];
          //upper[numUnderlyingSteps-1]            -= upper[numUnderlyingSteps-1] * ubc[3] / ubc[1];

          beagle::util::tridiagonalSolve( initialCondition, diag, upper, lower );
        }
      private:
        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
        beagle::real_2d_function_ptr_t m_Rate;
      };

      struct DupirePDESolver : public OneDimParabolicPDESolver
      {
        DupirePDESolver(const beagle::real_2d_function_ptr_t& convection,
                        const beagle::real_2d_function_ptr_t& diffusion) :
          m_Convection(convection),
          m_Diffusion(diffusion)
        { }
        virtual ~DupirePDESolver(void) = default;
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

            diag[j]  = 1. + 2. * diffusion * dTdXdX;
            upper[j] = - .5 * convection * dTdX - diffusion * dTdXdX;
            lower[j] =   .5 * convection * dTdX - diffusion * dTdXdX;
          }

          initialCondition[0] -= lower[0] * lbc[0];
          //diag[0]             -= lower[0] * lbc[2] / lbc[1];
          //upper[0]            -= lower[0] * lbc[3] / lbc[1];

          initialCondition[numUnderlyingSteps-1] -= upper[numUnderlyingSteps-1] * ubc[0];
          //diag[numUnderlyingSteps-1]             -= upper[numUnderlyingSteps-1] * ubc[2] / ubc[1];
          //upper[numUnderlyingSteps-1]            -= upper[numUnderlyingSteps-1] * ubc[3] / ubc[1];

          beagle::util::tridiagonalSolve( initialCondition, diag, upper, lower );
        }
      private:
        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
      };
    }


    beagle::parabolic_pde_solver_ptr_t
    OneDimParabolicPDESolver::formOneDimParabolicValuationPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                                    const beagle::real_2d_function_ptr_t& diffusion,
                                                                    const beagle::real_2d_function_ptr_t& rate,
                                                                    const beagle::real_2d_function_ptr_t& source)
    {
      return std::make_shared<impl::OneDimParabolicValuationPDESolver>(convection, diffusion, rate, source);
    }

    beagle::parabolic_pde_solver_ptr_t
    OneDimParabolicPDESolver::formOneDimFokkerPlanckPDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                              const beagle::real_2d_function_ptr_t& diffusion,
                                                              const beagle::real_2d_function_ptr_t& rate)
    {
      return std::make_shared<impl::OneDimFokkerPlanckPDESolver>(convection, diffusion, rate);
    }

    beagle::parabolic_pde_solver_ptr_t
    OneDimParabolicPDESolver::formDupirePDESolver(const beagle::real_2d_function_ptr_t& convection,
                                                  const beagle::real_2d_function_ptr_t& diffusion)
    {
      return std::make_shared<impl::DupirePDESolver>(convection, diffusion);
    }
  }
}
