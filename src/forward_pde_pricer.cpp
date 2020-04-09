#include "one_dim_pde_pricer.hpp"
#include "solver.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct OneDimensionalForwardPDEEuropeanOptionPricer : public OneDimensionalPDEOptionPricer,
                                                            public beagle::valuation::mixins::OptionValueCollectionProvider,
                                                            public beagle::valuation::mixins::CloneWithNewLocalVolatilitySurface
      {
        using two_dbl_t = std::pair<double, double>;

        OneDimensionalForwardPDEEuropeanOptionPricer( const FiniteDifferenceDetails& fdDetails,
                                                      const beagle::real_2d_function_ptr_t& volatility ) :
          OneDimensionalPDEOptionPricer( fdDetails ),
          m_Volatility( volatility )
        { }
        ~OneDimensionalForwardPDEEuropeanOptionPricer( void )
        { }
      public:
        virtual void formInitialOptionValueCollection( const beagle::payoff_ptr_t& payoff,
                                                       const beagle::dbl_vec_t& strikes,
                                                       beagle::dbl_vec_t& prices ) const
        {
          prices.resize( strikes.size() );
          std::transform( strikes.cbegin(),
                          strikes.cend(),
                          prices.begin(),
                          [&payoff, this](double strike) {return payoff->intrinsicValue(finiteDifferenceDetails().spot(), strike);}  );
        }
        virtual void optionValueCollection( double start,
                                            double end,
                                            const beagle::payoff_ptr_t& payoff,
                                            const beagle::dbl_vec_t& logStrikes,
                                            const beagle::dbl_vec_t& strikes,
                                            beagle::dbl_vec_t& prices ) const override
        {
          beagle::dbl_vec_t times;
          beagle::int_vec_t exDividendIndices;
          finiteDifferenceDetails().formTimeSteps( start, end, times, exDividendIndices );

          two_dbl_t boundaryStrikes = std::make_pair( std::exp(logStrikes.front()),
                                                      std::exp(logStrikes.back()) );

          int timeSteps = times.size();
          double deltaX = logStrikes[1] - logStrikes[0];

          int strikeSize = strikes.size();
          beagle::dbl_vec_t diag(strikeSize);
          beagle::dbl_vec_t lower(strikeSize);
          beagle::dbl_vec_t upper(strikeSize);

          auto it = finiteDifferenceDetails().dividends().cbegin();
          auto itEnd = finiteDifferenceDetails().dividends().cend();
          if (it != itEnd && it->first < start)
            ++it;

          auto jt = exDividendIndices.cbegin();
          auto jtEnd = exDividendIndices.cend();
          for (int i=0; i<timeSteps-1; ++i)
          {
            double thisTime = times[i+1];
            double deltaT = thisTime - times[i];
            for (int j=0; j<strikeSize; ++j)
            {
              double vol = m_Volatility->value(thisTime, strikes[j]);
              double volOverDeltaX = vol / deltaX;
              double volOverDeltaXSquared = volOverDeltaX * volOverDeltaX;
              double mu = finiteDifferenceDetails().rate() + .5 * vol * vol;
              double muOverDeltaX = mu / deltaX;
              diag[j]  = 1. + deltaT * volOverDeltaXSquared;
              upper[j] =   deltaT * .5 * (muOverDeltaX - volOverDeltaXSquared);
              lower[j] = - deltaT * .5 * (muOverDeltaX + volOverDeltaXSquared);
            }

            two_dbl_t boundaryValues = boundaryCondition( payoff,
                                                          boundaryStrikes,
                                                          end - thisTime );
            prices[0]            -= deltaT * lower[0] * boundaryValues.first;
            prices[strikeSize-1] -= deltaT * upper[strikeSize-1] * boundaryValues.second;

            beagle::util::tridiagonalSolve( prices, diag, upper, lower );

            /// Ex-dividend date
            if (jt != jtEnd && *jt == i)
            {
              beagle::real_function_ptr_t interpFunc = finiteDifferenceDetails().interpolation()->formFunction( strikes, prices );
              beagle::dbl_vec_t shiftedStrikes(strikes.size());

              double dividendAmount = it->second;
              std::transform( strikes.cbegin(),
                              strikes.cend(),
                              shiftedStrikes.begin(),
                              [this, dividendAmount](double strike) { 
                                return strike + finiteDifferenceDetails().dividendPolicy()->dividendAmount(finiteDifferenceDetails().spot(), dividendAmount);
                              } );

              std::transform( shiftedStrikes.cbegin(),
                              shiftedStrikes.cend(),
                              prices.begin(),
                              [&interpFunc](double strike) { 
                                return interpFunc->value(strike);
                              } );

              ++jt;
              ++it;
            }
          }
        }
        virtual double value( const beagle::product_ptr_t& product ) const override
        {
          if (auto pA = dynamic_cast<beagle::product::option::mixins::American*>(product.get()))
            throw(std::string("Cannot price an American option with forward PDE European option pricer!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();
          const beagle::payoff_ptr_t& payoff = pO->payoff();

          beagle::dbl_vec_t logStrikes;
          beagle::dbl_vec_t strikes;
          finiteDifferenceDetails().formStateVariableSteps( expiry, logStrikes, strikes );

          beagle::dbl_vec_t prices;
          formInitialOptionValueCollection( payoff, strikes, prices );
          optionValueCollection( 0., expiry, payoff, logStrikes, strikes, prices );
          // optionValueCollection( expiry/2., expiry, payoff, logStrikes, strikes, prices );

          beagle::real_function_ptr_t interpResult = finiteDifferenceDetails().interpolation()->formFunction( strikes, prices );
          return interpResult->value(strike);
        }
        virtual beagle::pricer_ptr_t createPricerWithNewLocalVolatilitySurface( const beagle::real_2d_function_ptr_t& vol ) const override
        {
          return std::make_shared<OneDimensionalForwardPDEEuropeanOptionPricer>( finiteDifferenceDetails(),
                                                                                 vol );
        }
      private:
        two_dbl_t boundaryCondition( const beagle::payoff_ptr_t& payoff,
                                     const two_dbl_t& boundaryStrikes,
                                     double timeToExpiry ) const
        {
          double minStrike = boundaryStrikes.first;
          double maxStrike = boundaryStrikes.second;

          double discounting = std::exp(-finiteDifferenceDetails().rate() * timeToExpiry);
          return std::make_pair( payoff->intrinsicValue(finiteDifferenceDetails().spot(), minStrike * discounting ),
                                 payoff->intrinsicValue(finiteDifferenceDetails().spot(), maxStrike * discounting ) );
        }
      private:
        beagle::real_2d_function_ptr_t m_Volatility;
      };

      struct OneDimForwardPDEEuroOptionPricer : public Pricer,
                                                public beagle::valuation::mixins::OneDimFokkerPlanck
      {
        OneDimForwardPDEEuroOptionPricer(const beagle::real_function_ptr_t& forward,
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
        virtual ~OneDimForwardPDEEuroOptionPricer( void )
        { }
      public:
        virtual double value(const beagle::product_ptr_t& product) const override
        {
          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          auto pA = dynamic_cast<beagle::product::option::mixins::American*>(product.get());

          double expiry = pO->expiry();
          double strike = pO->strike();
          const beagle::payoff_ptr_t& payoff = pO->payoff();

          // Set the initial condition and perform the forwardward induction
          beagle::dbl_vec_t stateVars;
          beagle::dbl_vec_t densities;
          formInitialCondition(expiry, stateVars, densities);
          evolve(0., expiry, stateVars, densities);
 
          double result(.0);
          int numStateVars = stateVars.size();
          double stateVarStep = (stateVars.back() - stateVars.front()) / (numStateVars - 1);
          double forward = m_Forward->value(expiry);
          for (int i=0; i<numStateVars; ++i)
          {
            result += stateVarStep * densities[i] * payoff->intrinsicValue(forward * std::exp(stateVars[i]), strike);
          }

          return m_Discounting->value(expiry) * result;
        }
        virtual void formInitialCondition(double expiry,
                                          beagle::dbl_vec_t& stateVars,
                                          beagle::dbl_vec_t& density) const override
        {
          double termForward = m_Forward->value(expiry);

          // Set up the state variable mesh
          int numStateVars = m_Settings.numberOfStateVariableSteps();
          if (numStateVars % 2 == 0)
            numStateVars += 1;

          double atmVol = m_Vol->value(expiry, termForward);
          int centralIndex = numStateVars / 2;
          double centralValue = 0.;
          double stateVarStep = m_Settings.numberOfGaussianStandardDeviations() * atmVol * std::sqrt(expiry) / centralIndex;

          stateVars.resize(numStateVars);
          for (int i=0; i<numStateVars; ++i)
            stateVars[i] = centralValue + (i - centralIndex) * stateVarStep;

          // Set the initial condition
          density.resize(numStateVars, 0.0);
          density[centralIndex] = 1. / stateVarStep;
        }
        virtual void evolve(double start,
                            double end,
                            const beagle::dbl_vec_t& stateVars,
                            beagle::dbl_vec_t& density) const override
        {
          beagle::real_2d_function_ptr_t convection = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [this](double time, double logMoneyness)
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
                            [this](double time, double logMoneyness)
                            {
                              double spot = m_Forward->value(time) * std::exp(logMoneyness);
                              return m_Rate->value(time, spot);
                            } );

          beagle::parabolic_pde_solver_ptr_t solver
            = beagle::math::OneDimParabolicPDESolver::formOneDimFokkerPlanckPDESolver(convection, diffusion, rate);

          // Perform the forward induction
          int numTimes = static_cast<int>((end - start) * m_Settings.numberOfTimeSteps());
          double timeStep = (end - start) / numTimes;
          for (int i=0; i<numTimes; ++i)
          {
            double end = start + (i+1) * timeStep;
            solver->evolve(end,
                           timeStep,
                           stateVars,
                           beagle::dbl_vec_t{0.0},
                           beagle::dbl_vec_t{0.0},
                           density);
          }
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
    Pricer::formOneDimensionalForwardPDEEuropeanOptionPricer( const FiniteDifferenceDetails& fdDetails,
                                                              const beagle::real_2d_function_ptr_t& volatility )
    {
      return std::make_shared<impl::OneDimensionalForwardPDEEuropeanOptionPricer>( fdDetails,
                                                                                   volatility );
    }

    beagle::pricer_ptr_t
    Pricer::formOneDimForwardPDEEuroOptionPricer(const beagle::real_function_ptr_t& forward,
                                                 const beagle::real_function_ptr_t& discounting,
                                                 const beagle::real_2d_function_ptr_t& drift,
                                                 const beagle::real_2d_function_ptr_t& volatility,
                                                 const beagle::real_2d_function_ptr_t& rate,
                                                 const beagle::valuation::OneDimFiniteDifferenceSettings& settings)
    {
      return std::make_shared<impl::OneDimForwardPDEEuroOptionPricer>( forward,
                                                                       discounting,
                                                                       drift,
                                                                       volatility,
                                                                       rate,
                                                                       settings );
    }
  }
}