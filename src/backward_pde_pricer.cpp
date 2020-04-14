#include "one_dim_pde_pricer.hpp"
#include "solver.hpp"
#include "bond.hpp"

#include <iterator>

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
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
          std::vector<beagle::two_dbl_t> cumAndExDivForwards;
          beagle::dividend_policy_ptr_t policy;
          auto pDS = dynamic_cast<beagle::math::mixins::ContainsDividends*>(forwardCurve().get());
          if (pDS)
          {
            auto schedule = pDS->dividendSchedule();
            auto it = std::lower_bound(schedule.cbegin(),
                                       schedule.cend(),
                                       expiry,
                                       [](const beagle::dividend_schedule_t::value_type& dividend,
                                          double value)
                                       { return std::get<0>(dividend) < value; });
            auto diff = std::distance(schedule.cbegin(), it);

            std::transform(schedule.cbegin(),
                           it,
                           std::back_inserter(dates),
                           [=](const beagle::dividend_schedule_t::value_type& item)
                           { return std::get<0>(item); });
            std::transform(schedule.cbegin(),
                           it,
                           std::back_inserter(dividends),
                           [=](const beagle::dividend_schedule_t::value_type& item)
                           { return std::make_pair(std::get<1>(item), std::get<2>(item)); });

            cumAndExDivForwards = pDS->cumAndExDividendForwards();
            cumAndExDivForwards.erase(cumAndExDivForwards.begin() + diff,
                                      cumAndExDivForwards.end());

            policy = pDS->dividendPolicy();

            if (cumAndExDivForwards.size() != dividends.size())
              throw(std::string(""));
            if (dates.size() - dividends.size() != 1U)
              throw(std::string(""));
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
          auto it = dates.crbegin();
          auto itEnd = dates.crend();
          auto jt = dividends.crbegin();
          auto jtEnd = dividends.crend();
          auto kt = cumAndExDivForwards.crbegin();
          double start = expiry;
          for ( ; it!=itEnd; ++it, ++jt, ++kt)
          {
            double end = *it;
            double interval = start - end;
            int numTimes = static_cast<int>(interval * finiteDifferenceSettings().numberOfTimeSteps());
            double timeStep = -interval / numTimes;
            for (int i=0; i<numTimes; ++i)
            {
              double thisTime = start - interval * (i+1) / numTimes;
              double df = discountCurve()->value(thisTime);
              double fwd = forwardCurve()->value(thisTime);
              if (i == numTimes-1 && jt != jtEnd)
                fwd = kt->second;

              double lbc = payoff->intrinsicValue( fwd * std::exp(stateVars.front() - stateVarStep), strike );
              double ubc = payoff->intrinsicValue( fwd * std::exp(stateVars.back()  + stateVarStep), strike );
              solver->evolve(thisTime,
                             timeStep,
                             stateVars,
                             beagle::dbl_vec_t{lbc},
                             beagle::dbl_vec_t{ubc},
                             prices);

              if (pA)
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

            // Take into account of dividends
            if (jt != jtEnd)
            {
              double cumForward = kt->first;
              double exForward = kt->second;

              beagle::dbl_vec_t spots(numStateVars);
              std::transform(moneynesses.cbegin(),
                             moneynesses.cend(),
                             spots.begin(),
                             [=](double moneyness)
                             { return exForward * moneyness; });

              beagle::real_function_ptr_t interpFunc
                = finiteDifferenceSettings().interpolationMethod()->formFunction( spots, prices );

              std::transform(moneynesses.cbegin(),
                             moneynesses.cend(),
                             spots.begin(),
                             [=](double moneyness)
                             { 
                               return policy->exDividendStockPrice(cumForward * moneyness * (1. - jt->first),
                                                                   jt->second);
                             } );

              std::transform( spots.cbegin(),
                              spots.cend(),
                              prices.begin(),
                              [&interpFunc](double spot)
                              { return interpFunc->value(spot); } );
            }

            start = end;
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