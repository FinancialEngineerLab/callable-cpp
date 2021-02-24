#include "pricer.hpp"
#include "util.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct ClosedFormSABREuropeanOptionPricer : public ClosedFormEuropeanOptionPricer
      {
        ClosedFormSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                           const beagle::real_function_ptr_t& discounting,
                                           const beagle::real_function_ptr_t& alpha,
                                           const beagle::real_function_ptr_t& beta,
                                           const beagle::real_function_ptr_t& rho,
                                           const beagle::real_function_ptr_t& nu) :
          ClosedFormEuropeanOptionPricer(forward, discounting),
          m_Alpha(alpha),
          m_Beta(beta),
          m_Rho(rho),
          m_Nu(nu)
        { }
      public:
        virtual double value( const beagle::product_ptr_t& product ) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();
          const auto& payoff = pO->payoff();
          double forward = forwardCurve()->value(expiry);
          double discounting = discountCurve()->value(expiry);

          return beagle::util::bsValue(strike, forward, expiry, impliedBlackScholesVolatility(product), payoff) * discounting;
        }
        virtual int numberOfParameters(void) const override
        {
          return 4;
        }
        virtual double impliedBlackScholesVolatility(const beagle::product_ptr_t& product) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();

          double forward = forwardCurve()->value(expiry);
          double alpha = m_Alpha->value(expiry);
          double beta = m_Beta->value(expiry);
          double rho = m_Rho->value(expiry);
          double nu = m_Nu->value(expiry);
          util::checkSABRParameters(alpha, beta, rho, nu, false);

          double result = 0.;

          if ( ! ( strike > 0.0 ) || ! ( forward > 0.0 ) )
            throw("Strike or forward in the SABR implied black volatility expansion must be positive");

          double halfOneMinusBeta = ( 1 - beta ) / 2.;
          double forwardTimesStrike = forward * strike;
          double forwardDividedByStrike = forward / strike;
          double tempOne = std::pow( forwardTimesStrike, halfOneMinusBeta );

          double factor = alpha / tempOne;
          double factorSquared = factor * factor;
          double halfOneMinusBetaSquared = halfOneMinusBeta * halfOneMinusBeta;
          double oneMinusBetaSquaredDividedByTwentyFour = halfOneMinusBetaSquared / 6.;
          result = factor * ( 1 + ( oneMinusBetaSquaredDividedByTwentyFour * factorSquared
                                  + .25 * beta * rho * nu * factor
                                  + ( 2. - 3. * rho * rho ) * nu * nu / 24. ) * expiry );

          if ( std::fabs( forward - strike ) > beagle::util::epsilon() )
          {
            double tempTwo = std::log( forwardDividedByStrike );
            double z = nu / alpha * tempOne * tempTwo;
            double chi = std::log( ( std::sqrt( 1. - 2. * rho * z + z * z ) + z - rho )
                                 / ( 1 - rho ) );

            result *= z / chi;

            double tempTwoSquared = tempTwo * tempTwo;
            double oneMinusBetaFourthDividedByOneNineTwoZero = halfOneMinusBetaSquared * halfOneMinusBetaSquared / 120.;
            double tempThree = 1. + oneMinusBetaSquaredDividedByTwentyFour * tempTwoSquared
                                  + oneMinusBetaFourthDividedByOneNineTwoZero * tempTwoSquared * tempTwoSquared;
            result /= tempThree;
          }

          return result;
        }
        beagle::real_function_ptr_coll_t modelParameters(void) const override
        {
          return {m_Alpha, m_Beta, m_Rho, m_Nu};
        }
        beagle::pricer_ptr_t updateModelParameters(const beagle::real_function_ptr_coll_t& params) const override
        {
          return Pricer::formClosedFormSABREuropeanOptionPricer(forwardCurve(),
                                                                discountCurve(),
                                                                params[0],
                                                                params[1],
                                                                params[2],
                                                                params[3]);
        }
      private:
        beagle::real_function_ptr_t m_Alpha;
        beagle::real_function_ptr_t m_Beta;
        beagle::real_function_ptr_t m_Rho;
        beagle::real_function_ptr_t m_Nu;
      };

      struct ClosedFormShiftedSABREuropeanOptionPricer : public ClosedFormEuropeanOptionPricer
      {
        ClosedFormShiftedSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                  const beagle::real_function_ptr_t& discounting,
                                                  const beagle::real_function_ptr_t& alpha,
                                                  const beagle::real_function_ptr_t& beta,
                                                  const beagle::real_function_ptr_t& rho,
                                                  const beagle::real_function_ptr_t& nu,
                                                  double shift) :
          ClosedFormEuropeanOptionPricer(forward, discounting),
          m_Alpha(alpha),
          m_Beta(beta),
          m_Rho(rho),
          m_Nu(nu),
          m_Shift(shift)
        { }
      public:
        virtual double value( const beagle::product_ptr_t& product ) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>(product.get());
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();
          const auto& payoff = pO->payoff();
          double forward = forwardCurve()->value(expiry);
          double discounting = discountCurve()->value(expiry);

          return beagle::util::bsValue(strike + m_Shift, forward + m_Shift, expiry, impliedBlackScholesVolatility(product), payoff) * discounting;
        }
        virtual int numberOfParameters(void) const override
        {
          return 4;
        }
        virtual double impliedBlackScholesVolatility(const beagle::product_ptr_t& product) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike() + m_Shift;

          double forward = forwardCurve()->value(expiry) + m_Shift;
          double alpha = m_Alpha->value(expiry);
          double beta = m_Beta->value(expiry);
          double rho = m_Rho->value(expiry);
          double nu = m_Nu->value(expiry);
          util::checkSABRParameters(alpha, beta, rho, nu, false);

          double result = 0.;

          if ( ! ( strike > 0.0 ) || ! ( forward > 0.0 ) )
            throw("Strike or forward in the SABR implied black volatility expansion must be positive");

          double halfOneMinusBeta = ( 1 - beta ) / 2.;
          double forwardTimesStrike = forward * strike;
          double forwardDividedByStrike = forward / strike;
          double tempOne = std::pow( forwardTimesStrike, halfOneMinusBeta );

          double factor = alpha / tempOne;
          double factorSquared = factor * factor;
          double halfOneMinusBetaSquared = halfOneMinusBeta * halfOneMinusBeta;
          double oneMinusBetaSquaredDividedByTwentyFour = halfOneMinusBetaSquared / 6.;
          result = factor * ( 1 + ( oneMinusBetaSquaredDividedByTwentyFour * factorSquared
                                  + .25 * beta * rho * nu * factor
                                  + ( 2. - 3. * rho * rho ) * nu * nu / 24. ) * expiry );

          if ( std::fabs( forward - strike ) > beagle::util::epsilon() )
          {
            double tempTwo = std::log( forwardDividedByStrike );
            double z = nu / alpha * tempOne * tempTwo;
            double chi = std::log( ( std::sqrt( 1. - 2. * rho * z + z * z ) + z - rho )
                                 / ( 1 - rho ) );

            result *= z / chi;

            double tempTwoSquared = tempTwo * tempTwo;
            double oneMinusBetaFourthDividedByOneNineTwoZero = halfOneMinusBetaSquared * halfOneMinusBetaSquared / 120.;
            double tempThree = 1. + oneMinusBetaSquaredDividedByTwentyFour * tempTwoSquared
                                  + oneMinusBetaFourthDividedByOneNineTwoZero * tempTwoSquared * tempTwoSquared;
            result /= tempThree;
          }

          return result;
        }
        beagle::real_function_ptr_coll_t modelParameters(void) const override
        {
          return {m_Alpha, m_Beta, m_Rho, m_Nu};
        }
        beagle::pricer_ptr_t updateModelParameters(const beagle::real_function_ptr_coll_t& params) const override
        {
          return Pricer::formClosedFormShiftedSABREuropeanOptionPricer(forwardCurve(),
                                                                       discountCurve(),
                                                                       params[0],
                                                                       params[1],
                                                                       params[2],
                                                                       params[3],
                                                                       m_Shift);
        }
      private:
        beagle::real_function_ptr_t m_Alpha;
        beagle::real_function_ptr_t m_Beta;
        beagle::real_function_ptr_t m_Rho;
        beagle::real_function_ptr_t m_Nu;
        double m_Shift;
      };
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                   const beagle::real_function_ptr_t& discounting,
                                                   const beagle::real_function_ptr_t& alpha,
                                                   const beagle::real_function_ptr_t& beta,
                                                   const beagle::real_function_ptr_t& rho,
                                                   const beagle::real_function_ptr_t& nu)
    {
      return std::make_shared<impl::ClosedFormSABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu);
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormShiftedSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                          const beagle::real_function_ptr_t& discounting,
                                                          const beagle::real_function_ptr_t& alpha,
                                                          const beagle::real_function_ptr_t& beta,
                                                          const beagle::real_function_ptr_t& rho,
                                                          const beagle::real_function_ptr_t& nu,
                                                          double shift)
    {
      return std::make_shared<impl::ClosedFormShiftedSABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, shift);
    }
  }
}
