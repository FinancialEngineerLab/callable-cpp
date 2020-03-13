#include "one_dim_pde_pricer.hpp"
#include "solver.hpp"

 #include <iostream>
// std::ofstream out("interpolation.txt");

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct OneDimensionalBackwardPDEOptionPricer : public OneDimensionalPDEOptionPricer,
                                                     public beagle::valuation::mixins::CloneWithNewLocalVolatilitySurface
      {
        using two_dbl_t = std::pair<double, double>;

        OneDimensionalBackwardPDEOptionPricer( const FiniteDifferenceDetails& fdDetails,
                                               const beagle::real_2d_function_ptr_t& volatility ) :
          OneDimensionalPDEOptionPricer( fdDetails ),
          m_Volatility(volatility)
        { }
        virtual ~OneDimensionalBackwardPDEOptionPricer( void )
        { }
      public:
        virtual double value( const beagle::product_ptr_t& product ) const override
        {
          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();
          const beagle::payoff_ptr_t& payoff = pO->payoff();

          beagle::dbl_vec_t times;
          beagle::int_vec_t exDividendIndices;
          beagle::dbl_vec_t logSpots;
          beagle::dbl_vec_t spots;
          formLatticeForBackwardValuation( expiry, times, exDividendIndices, logSpots, spots );

          auto pA = dynamic_cast<beagle::product::option::mixins::American*>(product.get());
          bool isAmerican( pA != nullptr );

          int spotSize = spots.size();

          // calculate terminal value for backward induction
          beagle::dbl_vec_t optionValues(spotSize);
          std::transform( spots.cbegin(),
                          spots.cend(),
                          optionValues.begin(),
                          [&payoff, strike](double spot) {return payoff->intrinsicValue(spot, strike);}  );

          two_dbl_t boundarySpots = std::make_pair( std::exp(logSpots.front()),
                                                    std::exp(logSpots.back()) );

          int timeSteps = times.size();
          double deltaX = logSpots[1] - logSpots[0];

          beagle::dbl_vec_t diag(spotSize);
          beagle::dbl_vec_t lower(spotSize);
          beagle::dbl_vec_t upper(spotSize);

          auto it = finiteDifferenceDetails().dividends().begin() + exDividendIndices.size() - 1;
          auto jt = exDividendIndices.crbegin();
          auto jtEnd = exDividendIndices.crend();
          for (int i=timeSteps-1; i>0; --i)
          {
            double thisTime = times[i-1];
            double deltaT = times[i] - thisTime;
            for (int j=0; j<spotSize; ++j)
            {
              double vol = m_Volatility->value(thisTime, spots[j]);
              double volOverDeltaX = vol / deltaX;
              double volOverDeltaXSquared = volOverDeltaX * volOverDeltaX;
              double mu = finiteDifferenceDetails().rate() - .5 * vol * vol;
              double muOverDeltaX = mu / deltaX;
              diag[j]  = 1. + deltaT * (volOverDeltaXSquared + finiteDifferenceDetails().rate());
              upper[j] = - deltaT * .5 * (muOverDeltaX + volOverDeltaXSquared);
              lower[j] =   deltaT * .5 * (muOverDeltaX - volOverDeltaXSquared);
            }

            two_dbl_t boundaryValues = boundaryCondition( payoff,
                                                          boundarySpots,
                                                          strike,
                                                          expiry - thisTime,
                                                          isAmerican );
            optionValues[0]          -= deltaT * lower[0] * boundaryValues.first;
            optionValues[spotSize-1] -= deltaT * upper[spotSize-1] * boundaryValues.second;

            beagle::util::tridiagonalSolve( optionValues, diag, upper, lower );

            if (isAmerican)
            {
              for (int j=0; j<spotSize; ++j)
              {
                optionValues[j] = std::max( payoff->intrinsicValue( spots[j], strike ), optionValues[j] );
              }
            }

            /// Ex-dividend date
            if (jt != jtEnd && *jt == i)
            {
              beagle::dbl_vec_t shiftedSpots(spots.cbegin(), spots.cend());
              beagle::real_function_ptr_t interpFunc = finiteDifferenceDetails().interpolation()->formFunction( spots, optionValues );

              // out << spots.size() << " " << optionValues.size() << std::endl;
              // for (int k=0; k<spots.size(); ++k)
              //   out  << spots[k] << " " << optionValues[k] << std::endl;

              double dividendAmount = it->second;
              std::transform( shiftedSpots.cbegin(),
                              shiftedSpots.cend(),
                              shiftedSpots.begin(),
                              [this, dividendAmount](double spot) { 
                                return finiteDifferenceDetails().dividendPolicy()->exDividendStockPrice(spot, dividendAmount);
                              } );

              std::transform( shiftedSpots.cbegin(),
                              shiftedSpots.cend(),
                              optionValues.begin(),
                              [&interpFunc](double spot) { 
                                return interpFunc->value(spot);
                              } );

              // out << std::endl << dividendAmount << std::endl;
              // for (int k=0; k<shiftedSpots.size(); ++k)
              //   out  << shiftedSpots[k] << " " << optionValues[k] << std::endl;
              // out << std::endl;

              ++jt;
              --it;
            }
          }

          beagle::real_function_ptr_t interpResult = finiteDifferenceDetails().interpolation()->formFunction( spots, optionValues );
          return interpResult->value(finiteDifferenceDetails().spot());
        }
        virtual beagle::pricer_ptr_t createPricerWithNewLocalVolatilitySurface( const beagle::real_2d_function_ptr_t& vol ) const override
        {
          return std::make_shared<OneDimensionalBackwardPDEOptionPricer>( finiteDifferenceDetails(),
                                                                          vol );
        }
      private:
        void formLatticeForBackwardValuation( double expiry,
                                              beagle::dbl_vec_t& times,
                                              beagle::int_vec_t& exDividendIndices,
                                              beagle::dbl_vec_t& logSpots,
                                              beagle::dbl_vec_t& spots ) const
        {
          finiteDifferenceDetails().formTimeSteps( 0., expiry, times, exDividendIndices );
          finiteDifferenceDetails().formStateVariableSteps( expiry, logSpots, spots );
        }
        two_dbl_t boundaryCondition( const beagle::payoff_ptr_t& payoff,
                                     const two_dbl_t& boundarySpots,
                                     double strike,
                                     double timeToExpiry,
                                     bool isAmerican ) const
        {
          double minSpot = boundarySpots.first;
          double maxSpot = boundarySpots.second;

          double adjustedStrike = strike * std::exp(-finiteDifferenceDetails().rate() * timeToExpiry);
          if (isAmerican && payoff->isPut())
            adjustedStrike = strike;

          return std::make_pair( payoff->intrinsicValue( minSpot, adjustedStrike ),
                                 payoff->intrinsicValue( maxSpot, adjustedStrike ) );
        }
      private:
        beagle::real_2d_function_ptr_t m_Volatility;
      };

      struct OneDimBackwardPDEOptionPricer : public Pricer
      {
        OneDimBackwardPDEOptionPricer(const beagle::real_function_ptr_t& forward,
                                      const beagle::real_function_ptr_t& discounting,
                                      const beagle::real_2d_function_ptr_t& drift,
                                      const beagle::real_2d_function_ptr_t& volatility,
                                      const beagle::real_2d_function_ptr_t& rate,
                                      const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
          m_Forward(forward),
          m_Discounting(discounting),
          m_Drift(drift),
          m_Vol(volatility),
          m_Rate(rate),
          m_Settings(settings)
        { }
        virtual ~OneDimBackwardPDEOptionPricer( void )
        { }
      public:
        virtual double value(const beagle::product_ptr_t& product) const
        {
          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          auto pA = dynamic_cast<beagle::product::option::mixins::American*>(product.get());

          double expiry = pO->expiry();
          double strike = pO->strike();
          const beagle::payoff_ptr_t& payoff = pO->payoff();

          // Calculate terminal forward and discount factor for later use
          double termForward = m_Forward->value(expiry);
          double termDF = m_Discounting->value(expiry);

          auto pDS = dynamic_cast<beagle::math::mixins::DividendSchedule*>(m_Forward.get());
          if (!pDS)
          {
            // Determine the convection and diffusion terms in the backward PDE
            beagle::real_function_ptr_t fundingRate = beagle::math::RealFunction::createUnaryFunction(
                              [this](double time)
                              {
                                double bump = 1e-6;
                                return std::log(m_Forward->value(time+bump) / m_Forward->value(time)) / bump;
                              } );
            beagle::real_function_ptr_t discountingRate = beagle::math::RealFunction::createUnaryFunction(
                              [this](double time)
                              {
                                double bump = 1e-6;
                                return - std::log(m_Discounting->value(time+bump) / m_Discounting->value(time)) / bump;
                              } );
            beagle::real_2d_function_ptr_t convection = beagle::math::RealTwoDimFunction::createBinaryFunction(
                              [this, fundingRate](double time, double logMoneyness)
                              {
                                double spot = m_Forward->value(time) * std::exp(logMoneyness);
                                double localVol = m_Vol->value(time, spot);
                                double drift = m_Drift->value(time, spot);
                                return drift - .5 * localVol * localVol;
                              } );
            beagle::real_2d_function_ptr_t diffusion = beagle::math::RealTwoDimFunction::createBinaryFunction(
                              [this](double time, double logMoneyness)
                              {
                                double spot = m_Forward->value(time) * std::exp(logMoneyness);
                                double localVol = m_Vol->value(time, spot);
                                return .5 * localVol * localVol;
                              } );
            beagle::real_2d_function_ptr_t rate = beagle::math::RealTwoDimFunction::createBinaryFunction(
                              [this, discountingRate](double time, double logMoneyness)
                              {
                                double spot = std::exp(logMoneyness);
                                return discountingRate->value(time) + m_Rate->value(time, spot);
                              } );

            beagle::parabolic_pde_solver_ptr_t solver
              = beagle::math::OneDimParabolicPDESolver::formOneDimParabolicValuationPDESolver(convection, diffusion,
                                                                                              beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.));

            // Set up the state variable mesh
            int numStateVars = m_Settings.numberOfStateVariableSteps();
            if (numStateVars % 2 == 0)
              numStateVars += 1;

            double atmVol = m_Vol->value(expiry, termForward);
            int centralIndex = numStateVars / 2;
            double centralValue = 0.;
            double stateVarStep = m_Settings.numberOfGaussianStandardDeviations() * atmVol * std::sqrt(expiry) / centralIndex;
            beagle::dbl_vec_t stateVars(numStateVars);
            for (int i=0; i<numStateVars; ++i)
              stateVars[i] = centralValue + (i - centralIndex) * stateVarStep;

            // Set the initial condition
            beagle::dbl_vec_t initialCondition(numStateVars);
            std::transform(stateVars.cbegin(),
                           stateVars.cend(),
                           initialCondition.begin(),
                           [=](double logMoneyness)
                           { return payoff->intrinsicValue(termForward * std::exp(logMoneyness), strike)
                                    * termDF; });

            // Perform the backward induction
            int numTimes = static_cast<int>(expiry * m_Settings.numberOfStateVariableSteps());
            double timeStep = (0. - expiry) / numTimes;
            for (int i=0; i<numTimes; ++i)
            {
              double end = expiry + (i+1) * timeStep / numTimes;
              double forward = m_Forward->value(end);
              double lbc = payoff->intrinsicValue( forward * std::exp(stateVars.front() - stateVarStep), strike );
              double ubc = payoff->intrinsicValue( forward * std::exp(stateVars.back()  + stateVarStep), strike );
              solver->evolve(end,
                             timeStep,
                             stateVars,
                             beagle::dbl_vec_t{lbc},
                             beagle::dbl_vec_t{ubc},
                             initialCondition);

              if (pA)
              {
                double df = m_Discounting->value(end);
                std::transform(initialCondition.begin(),
                               initialCondition.end(),
                               stateVars.cbegin(),
                               initialCondition.begin(),
                               [=](double continuation, double logSpot)
                               {
                                 return std::max(continuation,
                                                 payoff->intrinsicValue( forward * std::exp(logSpot), strike ) * termDF / df);
                               });
              }
            }

            return initialCondition[centralIndex];
          }
          else
            return 0.0;
        }
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        beagle::real_2d_function_ptr_t m_Drift;
        beagle::real_2d_function_ptr_t m_Vol;
        beagle::real_2d_function_ptr_t m_Rate;
        beagle::valuation::OneDimFiniteDifferenceSettings m_Settings;
      };
    }

    beagle::pricer_ptr_t
    Pricer::formOneDimensionalBackwardPDEOptionPricer( const FiniteDifferenceDetails& fdDetails,
                                                       const beagle::real_2d_function_ptr_t& diffusion)
    {
      return std::make_shared<impl::OneDimensionalBackwardPDEOptionPricer>( fdDetails,
                                                                            diffusion );
    }

    beagle::pricer_ptr_t
    Pricer::formOneDimBackwardPDEOptionPricer(const beagle::real_function_ptr_t& forward,
                                              const beagle::real_function_ptr_t& discounting,
                                              const beagle::real_2d_function_ptr_t& drift,
                                              const beagle::real_2d_function_ptr_t& volatility,
                                              const beagle::real_2d_function_ptr_t& rate,
                                              const beagle::valuation::OneDimFiniteDifferenceSettings& settings)
    {
      return std::make_shared<impl::OneDimBackwardPDEOptionPricer>( forward,
                                                                    discounting,
                                                                    drift,
                                                                    volatility,
                                                                    rate,
                                                                    settings );
    }
  }
}