#include "pricer.hpp"
#include "util.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "integration_method.hpp"

#include <stdio.h>

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      namespace
      {
        void checkCEVParameters(double beta,
                                double sigma,
                                bool isFreeBoundaryCEV)
        {
          double betaMax = isFreeBoundaryCEV ? .5 : 1.;
          const double& eps = util::epsilon();
          if (beta < eps || beta - betaMax > eps)
            throw("The beta parameter of the CEV model must be between zero and " + std::to_string(betaMax) + "!");
          if (sigma < eps)
            throw("The sigma parameter of the CEV model must be positive!");
        }
      }

      struct ClosedFormExactCEVEuropeanOptionPricer : public ClosedFormEuropeanOptionPricer
      {
        ClosedFormExactCEVEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                               const beagle::real_function_ptr_t& discounting,
                                               const beagle::real_function_ptr_t& beta,
                                               const beagle::real_function_ptr_t& sigma,
                                               const beagle::integration_method_ptr_t& quadMethod) :
          m_Forward(forward),
          m_Discounting(discounting),
          m_Beta(beta),
          m_Sigma(sigma),
          m_QuadMethod(quadMethod)
        { }
      public:
        double value(const beagle::product_ptr_t& product) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();
          double forward = m_Forward->value(expiry);
          double discouting = m_Discounting->value(expiry);
          double beta = m_Beta->value(expiry);
          double sigma = m_Sigma->value(expiry);
          checkCEVParameters(beta, sigma, false);
          const double& pi = util::pi();

          double result = 0.;
          double oneMinusBeta = 1 - beta;
          double nu = .5 / oneMinusBeta;
          double qK = std::pow( strike, oneMinusBeta ) / oneMinusBeta;
          double qF = std::pow( forward, oneMinusBeta ) / oneMinusBeta;
          double b = (qK * qK + qF * qF) / 2. / qK / qF;
          double tau = sigma * sigma * expiry;
          double factor = std::sqrt( forward * strike ) / pi;

          real_func_t fOne = [=](double theta) {
            double temp = b - std::cos(theta);
            return std::sin(nu * theta) * std::sin(theta) / temp
                   * std::exp( - qK * qF * temp / tau );
          };
          result += m_QuadMethod->quadrature(beagle::math::RealFunction::createUnaryFunction(fOne), 0, pi);

          real_func_t fTwo = [=](double theta) {
            double tanTheta = std::tan( theta );
            double temp = b + std::cosh(tanTheta);

            return std::sin(nu * pi) * std::exp(-nu * tanTheta) / temp
                    * std::exp(-qK * qF * temp / tau) * std::sinh(tanTheta)
                    * (1 + tanTheta * tanTheta);
          };
          result += m_QuadMethod->quadrature(beagle::math::RealFunction::createUnaryFunction(fTwo), 0, .499 * pi );

          if (pO->payoff()->isCall())
            result = factor * result + std::max(forward - strike, 0.);
          if (pO->payoff()->isPut())
            result = factor * result + std::max(strike - forward, 0.);

          return result * discouting;
        }
        int numberOfParameters(void) const override
        {
          return 2;
        }
        double impliedBlackScholesVolatility(const beagle::product_ptr_t& product) const override
        {
          return 0.0;
        }
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        beagle::real_function_ptr_t m_Beta;
        beagle::real_function_ptr_t m_Sigma;
        beagle::integration_method_ptr_t m_QuadMethod;
      };

      struct ClosedFormFreeBoundaryCEVEuropeanOptionPricer : public ClosedFormEuropeanOptionPricer
      {
        ClosedFormFreeBoundaryCEVEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                      const beagle::real_function_ptr_t& discounting,
                                                      const beagle::real_function_ptr_t& beta,
                                                      const beagle::real_function_ptr_t& sigma,
                                                      const beagle::integration_method_ptr_t& quadMethod) :
          m_Forward(forward),
          m_Discounting(discounting),
          m_Beta(beta),
          m_Sigma(sigma),
          m_QuadMethod(quadMethod)
        { }
      public:
        double value(const beagle::product_ptr_t& product) const override
        {
          auto pE = dynamic_cast<beagle::product::option::mixins::European*>(product.get());
          if (!pE)
            throw(std::string("Cannot value an option with non-European exercise style in closed form!"));

          auto pO = dynamic_cast<beagle::product::mixins::Option*>( product.get() );
          if (!pO)
            throw(std::string("The incoming product is not an option!"));

          double expiry = pO->expiry();
          double strike = pO->strike();
          double forward = m_Forward->value(expiry);
          double discouting = m_Discounting->value(expiry);
          double beta = m_Beta->value(expiry);
          double sigma = m_Sigma->value(expiry);
          checkCEVParameters(beta, sigma, true);
          const double& pi = util::pi();

          double result = 0.;
          double oneMinusBeta = 1 - beta;
          double eta = .5 / oneMinusBeta;
          double qK = std::pow( std::fabs(strike), oneMinusBeta ) / oneMinusBeta;
          double qF = std::pow( std::fabs(forward), oneMinusBeta ) / oneMinusBeta;
          double b = (qK * qK + qF * qF) / 2. / qK / qF;
          double tau = sigma * sigma * expiry;
          double factor = std::sqrt( std::fabs(forward * strike) ) / pi;

          real_func_t fOne = [=](double theta) {
            double temp = b - std::cos(theta);

            return std::sin(eta * theta) * std::sin(theta) / temp
                   * std::exp(-qK * qF * temp / tau);
          };
          if ( !(strike < util::epsilon()) )
            result += m_QuadMethod->quadrature(beagle::math::RealFunction::createUnaryFunction(fOne), 0, pi);

          real_func_t fTwo = [=](double theta) {
            double tanTheta = std::tan(theta);
            double temp = b + std::cosh(tanTheta);
            double common = std::sin(eta * pi) / temp
                             * (1 + tanTheta * tanTheta)
                             * std::exp( -qK * qF * temp / tau)
                             * std::sinh(tanTheta);
            if ( !(strike < util::epsilon()) )
              return common * std::cosh(eta * tanTheta);
            else
              return common * std::sinh(eta * tanTheta);
          };
          result += m_QuadMethod->quadrature(beagle::math::RealFunction::createUnaryFunction(fTwo), 0, .499 * pi);

          if (pO->payoff()->isCall())
            result = factor * result + std::max(forward - strike, 0.);
          if (pO->payoff()->isPut())
            result = factor * result + std::max(strike - forward, 0.);

          return result * discouting;
        }
        int numberOfParameters(void) const override
        {
          return 2;
        }
        double impliedBlackScholesVolatility(const beagle::product_ptr_t& product) const override
        {
          return 0.0;
        }
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        beagle::real_function_ptr_t m_Beta;
        beagle::real_function_ptr_t m_Sigma;
        beagle::integration_method_ptr_t m_QuadMethod;
      };
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormExactCEVEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                       const beagle::real_function_ptr_t& discounting,
                                                       const beagle::real_function_ptr_t& beta,
                                                       const beagle::real_function_ptr_t& sigma,
                                                       const beagle::integration_method_ptr_t& quadMethod)
    {
      return std::make_shared<impl::ClosedFormExactCEVEuropeanOptionPricer>(forward, discounting, beta, sigma, quadMethod);
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormFreeBoundaryCEVEuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                              const beagle::real_function_ptr_t& discounting,
                                                              const beagle::real_function_ptr_t& beta,
                                                              const beagle::real_function_ptr_t& sigma,
                                                              const beagle::integration_method_ptr_t& quadMethod)
    {
      return std::make_shared<impl::ClosedFormFreeBoundaryCEVEuropeanOptionPricer>(forward, discounting, beta, sigma, quadMethod);
    }
  }
}
