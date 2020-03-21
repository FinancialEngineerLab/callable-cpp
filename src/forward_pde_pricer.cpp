#include "one_dim_pde_pricer.hpp"
#include "solver.hpp"

#include <fstream>
std::ofstream out("distribution.txt");

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

          // for (auto price : prices)
          //   out << price << " ";
          // out << std::endl;

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

            // for (auto price : prices)
            //   out << price << " ";
            // out << std::endl;

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

      struct OneDimForwardPDEEuroOptionPricer : public Pricer
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
                                double spot =  std::exp(logMoneyness);
                                double localVol = m_Vol->value(time, spot);
                                double drift = m_Drift->value(time, spot);
                                return fundingRate->value(time) + drift - .5 * localVol * localVol;
                              } );
            beagle::real_2d_function_ptr_t diffusion = beagle::math::RealTwoDimFunction::createBinaryFunction(
                              [this](double time, double logMoneyness)
                              {
                                double spot = std::exp(logMoneyness);
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
              = beagle::math::OneDimParabolicPDESolver::formOneDimFokkerPlanckPDESolver(convection, diffusion, rate);

            // Set up the state variable mesh
            int numStateVars = m_Settings.numberOfStateVariableSteps();
            if (numStateVars % 2 == 0)
              numStateVars += 1;

            double atmVol = m_Vol->value(expiry, termForward);
            int centralIndex = numStateVars / 2;
            double centralValue = std::log(m_Forward->value(0.0));
            double stateVarStep = m_Settings.numberOfGaussianStandardDeviations() * atmVol * std::sqrt(expiry) / centralIndex;
            beagle::dbl_vec_t stateVars(numStateVars);
            for (int i=0; i<numStateVars; ++i)
              stateVars[i] = centralValue + (i - centralIndex) * stateVarStep;

            // Set the initial condition
            beagle::dbl_vec_t initialCondition(numStateVars, 0.0);
            initialCondition[centralIndex] = 1. / stateVarStep;

            // Perform the backward induction
            int numTimes = static_cast<int>(expiry * m_Settings.numberOfStateVariableSteps());
            double timeStep = expiry / numTimes;
            for (int i=0; i<numTimes; ++i)
            {
              double end = (i+1) * timeStep / numTimes;
              solver->evolve(end,
                             timeStep,
                             stateVars,
                             beagle::dbl_vec_t{0.0},
                             beagle::dbl_vec_t{0.0},
                             initialCondition);
            }
 
            double result(.0);
            for (int i=0; i<numStateVars; ++i)
            {
              result += stateVarStep * initialCondition[i] * payoff->intrinsicValue(std::exp(stateVars[i]), strike);
            }

            return result;
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