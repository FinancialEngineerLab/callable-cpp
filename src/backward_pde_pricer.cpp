#include "one_dim_pde_pricer.hpp"
#include "solver.hpp"
#include "bond.hpp"

#include <iostream>
#include <iterator>

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

      struct OneDimBackwardPDEOptionPricer : public OneDimParabolicPDEPricer
      {
        OneDimBackwardPDEOptionPricer(const beagle::real_function_ptr_t& forward,
                                      const beagle::real_function_ptr_t& discounting,
                                      const beagle::real_2d_function_ptr_t& drift,
                                      const beagle::real_2d_function_ptr_t& volatility,
                                      const beagle::real_2d_function_ptr_t& rate,
                                      const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
          OneDimParabolicPDEPricer(forward,
                                   discounting,
                                   drift,
                                   volatility,
                                   rate,
                                   beagle::math::RealTwoDimFunction::createTwoDimConstantFunction(0.0),
                                   settings)
        { }
        virtual ~OneDimBackwardPDEOptionPricer( void )
        { }
      public:
        virtual double value(const beagle::product_ptr_t& product) const
        {
          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();
          const beagle::payoff_ptr_t& payoff = pO->payoff();

          // Check American option
          auto pA = dynamic_cast<beagle::product::option::mixins::American*>(product.get());

          // Calculate terminal forward and discount factor for later use
          double termForward = forwardCurve()->value(expiry);
          double termDF = discountCurve()->value(expiry);

          // Check dividends and use ex-dividend dates in backward induction
          beagle::dbl_vec_t dates(1U, 0.0);
          std::vector<beagle::two_dbl_t> dividends;
          std::vector<beagle::two_dbl_t> cumAndExDivForward;
          beagle::dividend_policy_ptr_t policy;
          auto pDS = dynamic_cast<beagle::math::mixins::DividendSchedule*>(forwardCurve().get());
          if (pDS)
          {
            auto schedule = pDS->dividendSchedule();
            cumAndExDivForward = pDS->cumAndExDividendForwards();
            policy = pDS->dividendPolicy();

            dates.resize(schedule.size() + 1U);
            std::transform(schedule.cbegin(),
                           schedule.cend(),
                           dates.begin() + 1U,
                           [=](const beagle::dividend_schedule_t::value_type& item)
                           { return std::get<0>(item); });

            dividends.resize(schedule.size());
            std::transform(schedule.cbegin(),
                           schedule.cend(),
                           dividends.begin(),
                           [=](const beagle::dividend_schedule_t::value_type& item)
                           { return std::make_pair(std::get<1>(item), std::get<2>(item)); });
          }

          beagle::parabolic_pde_solver_ptr_t solver
            = beagle::math::OneDimParabolicPDESolver::formOneDimParabolicValuationPDESolver(
                                                      convectionCoefficient(),
                                                      diffusionCoefficient(),
                                                      rateCoefficient(),
                                                      sourceTerm());

          // Set up the state variable mesh
          int numStateVars = finiteDifferenceSettings().numberOfStateVariableSteps();
          if (numStateVars % 2 == 0)
            numStateVars += 1;

          double atmVol = volatilitySurface()->value(expiry, termForward);
          int centralIndex = numStateVars / 2;
          double centralValue = 0.;
          double stateVarStep = finiteDifferenceSettings().numberOfGaussianStandardDeviations() * atmVol * std::sqrt(expiry) / centralIndex;
          beagle::dbl_vec_t stateVars(numStateVars);
          for (int i=0; i<numStateVars; ++i)
            stateVars[i] = centralValue + (i - centralIndex) * stateVarStep;

          // Set the initial condition
          beagle::dbl_vec_t moneynesses(numStateVars);
          std::transform(stateVars.cbegin(),
                         stateVars.cend(),
                         moneynesses.begin(),
                         [=](double logMoneyness)
                         { return std::exp(logMoneyness); });

          beagle::dbl_vec_t prices(numStateVars);
          std::transform(moneynesses.cbegin(),
                         moneynesses.cend(),
                         prices.begin(),
                         [=](double moneyness)
                         { return termDF * payoff->intrinsicValue(termForward*moneyness, strike); });

          // Perform the backward induction
          int numTimes = static_cast<int>(expiry * finiteDifferenceSettings().numberOfTimeSteps());
          double timeStep = (0. - expiry) / numTimes;
          for (int i=0; i<numTimes; ++i)
          {
            double end = expiry * (numTimes - i - 1) / numTimes;
            double fwd = forwardCurve()->value(end);
            double df = discountCurve()->value(end);
            double lbc = payoff->intrinsicValue( fwd * std::exp(stateVars.front() - stateVarStep), strike );
            double ubc = payoff->intrinsicValue( fwd * std::exp(stateVars.back()  + stateVarStep), strike );
            solver->evolve(end,
                           timeStep,
                           stateVars,
                           beagle::dbl_vec_t{lbc},
                           beagle::dbl_vec_t{ubc},
                           prices);

            if (pA)
            {
              std::transform(prices.begin(),
                             prices.end(),
                             moneynesses.cbegin(),
                             prices.begin(),
                             [=](double continuation, double moneyness)
                             {
                               return std::max(continuation,
                                               df * payoff->intrinsicValue( fwd * moneyness, strike ));
                             });
            }
          }

          return prices[centralIndex];
        }
      };

      struct OneDimBackwardPDEBondPricer : public OneDimParabolicPDEPricer
      {
        OneDimBackwardPDEBondPricer(const beagle::real_function_ptr_t& forward,
                                    const beagle::real_function_ptr_t& discounting,
                                    const beagle::real_2d_function_ptr_t& drift,
                                    const beagle::real_2d_function_ptr_t& volatility,
                                    const beagle::real_2d_function_ptr_t& rate,
                                    const beagle::real_2d_function_ptr_t& recovery,
                                    const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
          OneDimParabolicPDEPricer(forward,
                                   discounting,
                                   drift,
                                   volatility,
                                   rate,
                                   recovery,
                                   settings)
        { }
        virtual ~OneDimBackwardPDEBondPricer( void )
        { }
      public:
        virtual double value(const beagle::product_ptr_t& product) const
        {
          // Extract bond cash flows
          auto pBond = dynamic_cast<beagle::product::mixins::Bond*>(product.get());
          if (!pBond)
            throw(std::string("The incoming product is not a bond!"));

          double notional = pBond->standardFaceValue();
          const beagle::coupon_flows_t& couponFlows = pBond->couponFlows();
          const beagle::notional_flows_t& notionalFlows = pBond->notionalFlows();

          if (notionalFlows.size() != 1U)
            throw(std::string("Cannot value a bond with notional structures!"));

          // Extract optionality features
          beagle::callable_schedule_t callSchedule;
          beagle::puttable_schedule_t putSchedule;
          beagle::real_function_ptr_t conversionRatio;
          bool isCallable(false);
          bool isPuttable(false);
          bool isConvertible(false);
          determineCallablePutableConvertibleFeatures(product,
                                                      callSchedule,
                                                      putSchedule,
                                                      conversionRatio,
                                                      isCallable,
                                                      isPuttable,
                                                      isConvertible);

          beagle::dbl_vec_t couponTimes(couponFlows.size() + 1U, 0.0);
          std::transform(couponFlows.cbegin(),
                         couponFlows.cend(),
                         couponTimes.begin() + 1,
                         [](const std::pair<double, double>& p)
                         { return p.first; });

          beagle::dbl_vec_t couponAmounts(couponFlows.size() + 1U, 0.0);
          std::transform(couponFlows.cbegin(),
                         couponFlows.cend(),
                         couponAmounts.begin() + 1,
                         [](const std::pair<double, double>& p)
                         { return p.second; });

          double expiry = notionalFlows.back().first;
          if (couponFlows.empty())
          {
            couponTimes.push_back(expiry);
            couponAmounts.push_back(0.0);
          }
          
          double termForward = forwardCurve()->value(expiry);
          double termDF = discountCurve()->value(expiry);

          beagle::parabolic_pde_solver_ptr_t solver
            = beagle::math::OneDimParabolicPDESolver::formOneDimParabolicValuationPDESolver(convectionCoefficient(),
                                                                                            diffusionCoefficient(),
                                                                                            rateCoefficient(),
                                                                                            sourceTerm());

          // Set up the state variable mesh
          int numStateVars = finiteDifferenceSettings().numberOfStateVariableSteps();
          if (numStateVars % 2 == 0)
            numStateVars += 1;

          double atmVol = volatilitySurface()->value(expiry, termForward);
          int centralIndex = numStateVars / 2;
          double centralValue = 0.;
          double stateVarStep = finiteDifferenceSettings().numberOfGaussianStandardDeviations() * atmVol * std::sqrt(expiry) / centralIndex;
          beagle::dbl_vec_t stateVars(numStateVars);
          for (int i=0; i<numStateVars; ++i)
            stateVars[i] = centralValue + (i - centralIndex) * stateVarStep;

          beagle::dbl_vec_t prices(numStateVars, termDF * notionalFlows.back().second);
          applyEarlyExerciseBoundaryConditions(stateVars,
                                               callSchedule,
                                               putSchedule,
                                               conversionRatio,
                                               isCallable,
                                               isPuttable,
                                               isConvertible,
                                               expiry,
                                               prices);

          int numCouponflows = couponTimes.size() - 1U;
          for (int j=0; j<numCouponflows; ++j)
          {
            double start = couponTimes[numCouponflows-j];
            double end = couponTimes[numCouponflows-j-1];

            // Pick up coupon
            double startDF = discountCurve()->value(start);
            std::transform(prices.begin(),
                           prices.end(),
                           prices.begin(),
                           [=](double price)
                           { return price + startDF * couponAmounts[numCouponflows-j]; });
            
            // Perform the backward induction
            int numTimes = static_cast<int>((start - end) * finiteDifferenceSettings().numberOfTimeSteps());
            double timeStep = (end - start) / numTimes;
            for (int i=0; i<numTimes; ++i)
            {
              double thisTime = start + (end - start) * (i+1) / numTimes;
              double thisDF = discountCurve()->value(thisTime);
              double lbc = thisDF * notional;
              double ubc = thisDF * notional;
              solver->evolve(thisTime,
                             timeStep,
                             stateVars,
                             beagle::dbl_vec_t{lbc},
                             beagle::dbl_vec_t{ubc},
                             prices);

              applyEarlyExerciseBoundaryConditions(stateVars,
                                                   callSchedule,
                                                   putSchedule,
                                                   conversionRatio,
                                                   isCallable,
                                                   isPuttable,
                                                   isConvertible,
                                                   thisTime,
                                                   prices);
            }
          }

          return prices[centralIndex];
        }
      private:
        void determineCallablePutableConvertibleFeatures( const beagle::product_ptr_t& product,
                                                          beagle::callable_schedule_t& callSchedule,
                                                          beagle::puttable_schedule_t& putSchedule,
                                                          beagle::real_function_ptr_t& convRatio,
                                                          bool& isCallable,
                                                          bool& isPuttable,
                                                          bool& isConvertible ) const
        {
          auto pCall = dynamic_cast<beagle::product::bond::mixins::Callable*>(product.get());
          if (pCall)
          {
            callSchedule = pCall->callSchedule();
            isCallable = true;

            if (callSchedule.size() > 1U)
              throw(std::string("Multiple call period is not supported currently!"));
          }
          
          auto pPut = dynamic_cast<beagle::product::bond::mixins::Puttable*>(product.get());
          if (pPut)
          {
            putSchedule = pPut->putSchedule();
            isPuttable = true;
          }
          
          auto pConv = dynamic_cast<beagle::product::bond::mixins::Convertible*>(product.get());
          if (pConv)
          {
            convRatio = pConv->conversionRatio();
            isConvertible = true;
          }
        }
        void applyEarlyExerciseBoundaryConditions(const beagle::dbl_vec_t& stateVars,
                                                  const beagle::callable_schedule_t& callSchedule,
                                                  const beagle::puttable_schedule_t& putSchedule,
                                                  const beagle::real_function_ptr_t& convRatio,
                                                  bool isCallable,
                                                  bool isPuttable,
                                                  bool isConvertible,
                                                  double time,
                                                  beagle::dbl_vec_t& prices) const
        {
          int sz = prices.size();
          double forward = forwardCurve()->value(time);
          double df = discountCurve()->value(time);

          // Put date
          if (isPuttable)
          {
            auto jt = std::find_if(putSchedule.cbegin(),
                                   putSchedule.cend(),
                                   [time](const puttable_schedule_t::value_type& val)
                                   { return std::fabs(val.first - time) < .001; });

            if (jt != putSchedule.cend())
            {
              for (int k=0; k<sz; ++k)
                prices[k] = std::max(prices[k], df * jt->second);
            }
          }

          if (isCallable)
          {
            double callPeriodStart = std::get<0>(callSchedule.front());
            double callPeriodEnd   = std::get<2>(callSchedule.front());
            const beagle::real_function_ptr_t& callPrice = std::get<1>(callSchedule.front());

            // In the callable period
            if (time - callPeriodStart >= .0 && time - callPeriodEnd <= .0)
            {
              // Callable and convertible
              if (isConvertible)
              {
                for (int k=0; k<sz; ++k)
                {
                  double conversionValue = convRatio->value(time) * forward * std::exp(stateVars[k]);
                  double callValue = callPrice->value(time);
                  prices[k] = std::min(prices[k],
                                       df * std::max(conversionValue, callValue));
                }
              }
              // Callable only
              else
              {
                for (int k=0; k<sz; ++k)
                {
                  double callValue = callPrice->value(time);
                  prices[k] = std::min(prices[k], df * callValue);
                }
              }
            }
          }
          
          if (isConvertible)
          {
            for (int k=0; k<sz; ++k)
            {
              double conversionValue = convRatio->value(time) * forward * std::exp(stateVars[k]);
              prices[k] = std::max(prices[k], df * conversionValue);
            }
          }
        }
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
    
    beagle::pricer_ptr_t
    Pricer::formOneDimBackwardPDEBondPricer(const beagle::real_function_ptr_t& forward,
                                            const beagle::real_function_ptr_t& discounting,
                                            const beagle::real_2d_function_ptr_t& drift,
                                            const beagle::real_2d_function_ptr_t& volatility,
                                            const beagle::real_2d_function_ptr_t& rate,
                                            const beagle::real_2d_function_ptr_t& recovery,
                                            const beagle::valuation::OneDimFiniteDifferenceSettings& settings)
    {
      return std::make_shared<impl::OneDimBackwardPDEBondPricer>( forward,
                                                                  discounting,
                                                                  drift,
                                                                  volatility,
                                                                  rate,
                                                                  recovery,
                                                                  settings );
    }
  }
}