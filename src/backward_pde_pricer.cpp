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
        using two_dbl_t = std::pair<beagle::dbl_t, beagle::dbl_t>;

        OneDimensionalBackwardPDEOptionPricer( const FiniteDifferenceDetails& fdDetails,
                                               const beagle::real_2d_function_ptr_t& volatility ) :
          OneDimensionalPDEOptionPricer( fdDetails ),
          m_Volatility(volatility)
        { }
        virtual ~OneDimensionalBackwardPDEOptionPricer( void )
        { }
      public:
        virtual beagle::dbl_t value( const beagle::product_ptr_t& product ) const override
        {
          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          beagle::dbl_t expiry = pO->expiry();
          beagle::dbl_t strike = pO->strike();
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
                          [&payoff, strike](beagle::dbl_t spot) {return payoff->intrinsicValue(spot, strike);}  );

          two_dbl_t boundarySpots = std::make_pair( std::exp(logSpots.front()),
                                                    std::exp(logSpots.back()) );

          int timeSteps = times.size();
          beagle::dbl_t deltaX = logSpots[1] - logSpots[0];

          beagle::dbl_vec_t diag(spotSize);
          beagle::dbl_vec_t lower(spotSize);
          beagle::dbl_vec_t upper(spotSize);

          auto it = finiteDifferenceDetails().dividends().begin() + exDividendIndices.size() - 1;
          auto jt = exDividendIndices.crbegin();
          auto jtEnd = exDividendIndices.crend();
          for (int i=timeSteps-1; i>0; --i)
          {
            beagle::dbl_t thisTime = times[i-1];
            beagle::dbl_t deltaT = times[i] - thisTime;
            for (int j=0; j<spotSize; ++j)
            {
              beagle::dbl_t vol = m_Volatility->value(thisTime, spots[j]);
              beagle::dbl_t volOverDeltaX = vol / deltaX;
              beagle::dbl_t volOverDeltaXSquared = volOverDeltaX * volOverDeltaX;
              beagle::dbl_t mu = finiteDifferenceDetails().rate() - .5 * vol * vol;
              beagle::dbl_t muOverDeltaX = mu / deltaX;
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

              beagle::dbl_t dividendAmount = it->second;
              std::transform( shiftedSpots.cbegin(),
                              shiftedSpots.cend(),
                              shiftedSpots.begin(),
                              [this, dividendAmount](beagle::dbl_t spot) { 
                                return finiteDifferenceDetails().dividendPolicy()->exDividendStockPrice(spot, dividendAmount);
                              } );

              std::transform( shiftedSpots.cbegin(),
                              shiftedSpots.cend(),
                              optionValues.begin(),
                              [&interpFunc](beagle::dbl_t spot) { 
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
        void formLatticeForBackwardValuation( beagle::dbl_t expiry,
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
                                     beagle::dbl_t strike,
                                     beagle::dbl_t timeToExpiry,
                                     bool isAmerican ) const
        {
          beagle::dbl_t minSpot = boundarySpots.first;
          beagle::dbl_t maxSpot = boundarySpots.second;

          beagle::dbl_t adjustedStrike = strike * std::exp(-finiteDifferenceDetails().rate() * timeToExpiry);
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
        virtual beagle::dbl_t value(const beagle::product_ptr_t& product) const
        {
          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          auto pA = dynamic_cast<beagle::product::option::mixins::American*>(product.get());

          beagle::dbl_t expiry = pO->expiry();
          beagle::dbl_t strike = pO->strike();
          const beagle::payoff_ptr_t& payoff = pO->payoff();

          // Calculate terminal forward and discount factor for later use
          beagle::dbl_t termForward = m_Forward->value(expiry);
          beagle::dbl_t termDF = m_Discounting->value(expiry);

          auto pDS = dynamic_cast<beagle::math::mixins::DividendSchedule*>(m_Forward.get());
          if (!pDS)
          {
            // Determine the convection and diffusion terms in the backward PDE
            beagle::real_function_ptr_t fundingRate = beagle::math::RealFunction::createUnaryFunction(
                              [this](beagle::dbl_t time)
                              {
                                beagle::dbl_t bump = 1e-6;
                                return std::log(m_Forward->value(time+bump) / m_Forward->value(time)) / bump;
                              } );
            beagle::real_function_ptr_t discountingRate = beagle::math::RealFunction::createUnaryFunction(
                              [this](beagle::dbl_t time)
                              {
                                beagle::dbl_t bump = 1e-6;
                                return - std::log(m_Discounting->value(time+bump) / m_Discounting->value(time)) / bump;
                              } );
            beagle::real_2d_function_ptr_t convection = beagle::math::RealTwoDimFunction::createBinaryFunction(
                              [this, fundingRate](beagle::dbl_t time, beagle::dbl_t logMoneyness)
                              {
                                beagle::dbl_t spot =  std::exp(logMoneyness);
                                beagle::dbl_t localVol = m_Vol->value(time, spot);
                                beagle::dbl_t drift = m_Drift->value(time, spot);
                                return fundingRate->value(time) + drift - .5 * localVol * localVol;
                              } );
            beagle::real_2d_function_ptr_t diffusion = beagle::math::RealTwoDimFunction::createBinaryFunction(
                              [this](beagle::dbl_t time, beagle::dbl_t logMoneyness)
                              {
                                beagle::dbl_t spot = std::exp(logMoneyness);
                                beagle::dbl_t localVol = m_Vol->value(time, spot);
                                return .5 * localVol * localVol;
                              } );
            beagle::real_2d_function_ptr_t rate = beagle::math::RealTwoDimFunction::createBinaryFunction(
                              [this, discountingRate](beagle::dbl_t time, beagle::dbl_t logMoneyness)
                              {
                                beagle::dbl_t spot = std::exp(logMoneyness);
                                return discountingRate->value(time) + m_Rate->value(time, spot);
                              } );

            beagle::parabolic_pde_solver_ptr_t solver
              = beagle::math::OneDimParabolicPDESolver::formOneDimParabolicValuationPDESolver(convection, diffusion, m_Rate);

            // Set up the state variable mesh
            int numStateVars = m_Settings.numberOfStateVariableSteps();
            if (numStateVars % 2 == 0)
              numStateVars += 1;

            beagle::dbl_t atmVol = m_Vol->value(expiry, termForward);
            int centralIndex = numStateVars / 2;
            beagle::dbl_t centralValue = std::log(m_Forward->value(0.0));
            beagle::dbl_t stateVarStep = m_Settings.numberOfGaussianStandardDeviations() * atmVol * std::sqrt(expiry) / centralIndex;
            beagle::dbl_vec_t stateVars(numStateVars);
            for (int i=0; i<numStateVars; ++i)
              stateVars[i] = centralValue + (i - centralIndex) * stateVarStep;

            // Set the initial condition
            beagle::dbl_vec_t initialCondition(numStateVars);
            std::transform(stateVars.cbegin(),
                           stateVars.cend(),
                           initialCondition.begin(),
                           [=](beagle::dbl_t logMoneyness)
                           { return payoff->intrinsicValue(std::exp(logMoneyness), strike)
                                    * termDF; });

            // Perform the backward induction
            int numTimes = static_cast<int>(expiry * m_Settings.numberOfStateVariableSteps());
            beagle::dbl_t timeStep = (0. - expiry) / numTimes;
            for (int i=0; i<numTimes; ++i)
            {
              beagle::dbl_t end = expiry + (i+1) * timeStep / numTimes;
              beagle::dbl_t lbc = payoff->intrinsicValue( std::exp(stateVars.front() - stateVarStep), strike );
              beagle::dbl_t ubc = payoff->intrinsicValue( std::exp(stateVars.back() + stateVarStep), strike );
              solver->evolve(end,
                             timeStep,
                             stateVars,
                             beagle::dbl_vec_t{lbc},
                             beagle::dbl_vec_t{ubc},
                             initialCondition);

              if (pA)
              {
                beagle::dbl_t df = m_Discounting->value(end);
                std::transform(initialCondition.begin(),
                               initialCondition.end(),
                               stateVars.cbegin(),
                               initialCondition.begin(),
                               [=](beagle::dbl_t continuation, beagle::dbl_t logSpot)
                               {
                                 return std::max(continuation,
                                                 payoff->intrinsicValue( std::exp(logSpot), strike ) * df);
                               });
              }
            }

            //for (int i=0; i<numStateVars; ++i)
            //  std::cout << "\n" << stateVars[i] << ":\t\t" << initialCondition[i];
            //std::cout << "\n\n";

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