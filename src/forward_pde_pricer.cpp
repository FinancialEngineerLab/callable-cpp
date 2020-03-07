#include "one_dim_pde_pricer.hpp"

// #include <fstream>
// std::ofstream out("interpolation.txt");

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

      struct OneDimensionalForwardPDEPricer : public Pricer
      {
        OneDimensionalForwardPDEPricer(double spot,
                                       const beagle::real_2d_function_ptr_t& drift,
                                       const beagle::real_2d_function_ptr_t& volatility,
                                       const beagle::real_function_ptr_t& rate,
                                       const beagle::discrete_dividend_schedule_t& dividends,
                                       const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
          m_Spot(spot),
          m_Drift(drift),
          m_Vol(volatility),
          m_Rate(rate),
          m_Dividends(dividends),
          m_Settings(settings)
        { }
        virtual ~OneDimensionalForwardPDEPricer( void )
        { }
      public:
        virtual double value(const beagle::product_ptr_t& product) const
        {
          return 0.0;
        }
      private:
        double m_Spot;
        beagle::real_2d_function_ptr_t m_Drift;
        beagle::real_2d_function_ptr_t m_Vol;
        beagle::real_function_ptr_t m_Rate;
        beagle::discrete_dividend_schedule_t m_Dividends;
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
    Pricer::formOneDimensionalForwardPDEPricer(double spot,
                                               const beagle::real_2d_function_ptr_t& drift,
                                               const beagle::real_2d_function_ptr_t& volatility,
                                               const beagle::real_function_ptr_t& rate,
                                               const beagle::discrete_dividend_schedule_t& dividends,
                                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings)
    {
      return std::make_shared<impl::OneDimensionalForwardPDEPricer>( spot,
                                                                     drift,
                                                                     volatility,
                                                                     rate,
                                                                     dividends,
                                                                     settings );
    }
  }
}