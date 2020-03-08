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
        virtual beagle::dbl_vec_t lowerBoundaryCoefficients( void ) const
        {
          return beagle::dbl_vec_t{0., 1., 0., 0.};
        }
        virtual beagle::dbl_vec_t upperBoundaryCoefficients( void ) const
        {
          return beagle::dbl_vec_t{0., 0., 0., 1.};
        }
      };
      
      struct LinearBoundaryCondition : public BoundaryCondition
      {
        virtual ~LinearBoundaryCondition( void )
        { }
      public:
        virtual beagle::dbl_vec_t lowerBoundaryCoefficients( void ) const
        {
          return beagle::dbl_vec_t{0., 1., -2., 1.};
        }
        virtual beagle::dbl_vec_t upperBoundaryCoefficients( void ) const
        {
          return beagle::dbl_vec_t{0., 1., -2., 1.};
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
                                                       const beagle::real_2d_function_ptr_t& convection,
                                                       const beagle::real_2d_function_ptr_t& diffusion,
                                                       const beagle::real_function_ptr_t& rate,
                                                       const beagle::discrete_dividend_schedule_t& dividends,
                                                       const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
      m_Spot(spot),
      m_Convection(convection),
      m_Diffusion(diffusion),
      m_Rate(rate),
      m_Dividends(dividends),
      m_Settings(settings)
    { }

    beagle::real_function_ptr_t
    OneDimParabolicPDESolver::evolve(double start,
                                     double end,
                                     const beagle::boundary_condition_ptr_t& lowerBoundaryCondition,
                                     const beagle::boundary_condition_ptr_t& upperBoundaryCondition,
                                     const beagle::real_function_ptr_t& initialCondition) const
    {
      // Form time steps
      beagle::dbl_vec_t times;
      double expiry = end - start;
      int numTimeSteps = static_cast<int>(std::floor(expiry * m_Settings.numberOfTimeSteps()));
      times.resize(numTimeSteps + 1);
      for (int i=0; i<numTimeSteps+1; ++i)
        times[i] = start + i * expiry / numTimeSteps;

      // Form underlying steps
      // Notice we are writing the PDE in the log price space, the drift and the 
      // volatility should be evaluated at the exponential of the lattice points

      // Two extra points are reserved for boundary conditions
      int numUnderlyingSteps = m_Settings.numberOfStateVariableSteps();
      beagle::dbl_vec_t underlyings(numUnderlyingSteps);

      double center = std::log(m_Spot) + expiry * m_Convection->value(end, m_Spot);
      double spread = m_Settings.numberOfGaussianStandardDeviations()
                    * std::sqrt(2. * m_Diffusion->value(end, m_Spot) * end);

      int centralIndex = numUnderlyingSteps / 2;
      double deltaX = spread / centralIndex;
      for (int i=0; i<numUnderlyingSteps; ++i)
        underlyings[i] = center + (i - centralIndex) * deltaX;

      // Discretize the intial condition
      beagle::dbl_vec_t prices(numUnderlyingSteps);
      beagle::dbl_vec_t svs(numUnderlyingSteps);
      std::transform(underlyings.cbegin(),
                     underlyings.cend(),
                     svs.begin(),
                     [](double z)
                     { return std::exp(z); });
      std::transform(svs.cbegin(),
                     svs.cend(),
                     prices.begin(),
                     [&initialCondition](double z)
                     { return initialCondition->value(z); });

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
          double convection = m_Convection->value(thisTime, svs[j]);
          double diffusion = m_Diffusion->value(thisTime, svs[j]);
          double rate = m_Rate->value(thisTime);

          diag[j]  = 1. - rate * deltaT - 2. * diffusion * dTdXdX;
          upper[j] =   .5 * convection * dTdX + diffusion * dTdXdX;
          lower[j] = - .5 * convection * dTdX + diffusion * dTdXdX;
        }

        beagle::dbl_vec_t lbc = lowerBoundaryCondition->lowerBoundaryCoefficients();
        prices[0] -= lower[0] * lbc[0] / lbc[1];
        diag[0]   -= lower[0] * lbc[2] / lbc[1];
        upper[0]  -= lower[0] * lbc[3] / lbc[1];

        beagle::dbl_vec_t ubc = upperBoundaryCondition->upperBoundaryCoefficients();
        prices[numUnderlyingSteps-1] -= upper[numUnderlyingSteps-1] * ubc[0] / ubc[3];
        diag[numUnderlyingSteps-1]   -= upper[numUnderlyingSteps-1] * ubc[1] / ubc[3];
        upper[numUnderlyingSteps-1]  -= upper[numUnderlyingSteps-1] * ubc[2] / ubc[3];

        // for (auto price : prices)
        //   out << price << " ";
        // out << std::endl;

        beagle::util::tridiagonalSolve( prices, diag, upper, lower );
      }

      return m_Settings.interpolationMethod()->formFunction(svs, prices);
    }
  }
}