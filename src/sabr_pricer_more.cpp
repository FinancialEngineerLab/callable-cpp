#include "pricer.hpp"
#include "util.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "real_function.hpp"
#include "integration_method.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct ClosedFormMappedZeroCorrelationSABREuropeanOptionPricer : public ClosedFormEuropeanOptionPricer,
                                                                       public mixins::MappedZeroCorrelationSABRParameters
      {
        ClosedFormMappedZeroCorrelationSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                const beagle::real_function_ptr_t& discounting,
                                                                const beagle::real_function_ptr_t& alpha,
                                                                const beagle::real_function_ptr_t& beta,
                                                                const beagle::real_function_ptr_t& rho,
                                                                const beagle::real_function_ptr_t& nu,
                                                                bool useApproximateKernel,
                                                                const beagle::integration_method_ptr_t& quadMethod) :
          m_Alpha(alpha),
          m_Beta(beta),
          m_Rho(rho),
          m_Nu(nu),
          m_KernelMethod(useApproximateKernel),
          m_QuadMethod(quadMethod)
        { }
      public:
        int numberOfParameters( void ) const override
        {
          return 4;
        }
        beagle::real_function_ptr_coll_t modelParameters( void ) const override
        {
          return {m_Alpha, m_Beta, m_Rho, m_Nu};
        }
        double impliedBlackScholesVolatility(const beagle::product_ptr_t& product) const override
        {
          return 0.0;
        }
      public:
        double effectiveZeroCorrelationAlpha(double strike,
                                             double forward,
                                             double expiry)  const override
        {
          double alpha = m_Alpha->value(expiry);
          double beta = m_Beta->value(expiry);
          double rho = m_Rho->value(expiry);
          double nu = m_Nu->value(expiry);
          const double& eps = beagle::util::epsilon();
          const double& pi = beagle::util::pi();

          if (std::fabs( rho ) < eps)
            return alpha;
          else
          {
            double mappedAlphaZerothOrder;
            if ( std::fabs( forward - strike ) > eps )
            {
              double mappedBeta = effectiveZeroCorrelationBeta(strike, forward, expiry);
              double mappedNu = effectiveZeroCorrelationNu(strike, forward, expiry);

              double effectiveStrike = std::fabs(doEffectiveStrike( strike, forward ));
              double effectiveForward = std::fabs( forward );

              double oneMinusBeta = 1 - beta;
              double qK = std::pow( effectiveStrike, oneMinusBeta ) / oneMinusBeta;
              double qF = std::pow( effectiveForward, oneMinusBeta ) / oneMinusBeta;
              double deltaQ = qK - qF;

              double oneMinusMappedBeta = 1. - mappedBeta;
              double qKMapped = std::pow( effectiveStrike, oneMinusMappedBeta ) / oneMinusMappedBeta;
              double qFMapped = std::pow( effectiveForward, oneMinusMappedBeta ) / oneMinusMappedBeta;
              double deltaQMapped = qKMapped - qFMapped;

              double alphaMinSquared = nu*nu*deltaQ*deltaQ + 2*rho*nu*deltaQ*alpha + alpha*alpha;
              double alphaMin = std::sqrt( alphaMinSquared );
              double arg = (alphaMin + rho*alpha + nu*deltaQ) / (1.+rho) / alpha;
              double phi = std::pow( arg, mappedNu/nu );

              if (mappedNu < eps)
                mappedAlphaZerothOrder = nu * deltaQMapped / std::log( arg );
              else
                mappedAlphaZerothOrder = 2.*phi*deltaQMapped*mappedNu / (phi*phi - 1);

              double alphaMinMapped =
                    std::sqrt(mappedNu*mappedNu*deltaQMapped*deltaQMapped + mappedAlphaZerothOrder*mappedAlphaZerothOrder);
              double rootOneMinusRhoSquared = std::sqrt(1. - rho*rho);
              double L = alphaMin / qK / nu / rootOneMinusRhoSquared;
              double uNaught = (deltaQ * nu * rho + alpha - alphaMin) / deltaQ / nu / rootOneMinusRhoSquared;

              double I = 2. * uNaught / (1 + uNaught);
              if (L - 1. < 0)
              {
                double rootOneMinusLSquared = std::sqrt(1. - L * L);
                double tempOne = (uNaught + L) / rootOneMinusLSquared;
                double tempTwo = L / rootOneMinusLSquared;
                I = 2. * (std::atan(tempOne) - std::atan(tempTwo)) / rootOneMinusLSquared;
              }
              else if (L - 1. > 0)
              {
                double rootLSquaredMinusOne = std::sqrt(L * L - 1.);
                double tempOne = L + rootLSquaredMinusOne;
                double tempTwo = L - rootLSquaredMinusOne;
                I = (std::log(1. + uNaught * tempOne) - std::log(1. + uNaught * tempTwo)) / rootLSquaredMinusOne;
              }

              double B = -.5 * beta / oneMinusBeta * rho / rootOneMinusRhoSquared
                       * (pi - std::acos( -(deltaQ * nu + alpha * rho) / alphaMin ) - std::acos( rho ) - I);

              double mappedAlphaFirstOrder = mappedNu * mappedNu / std::log(phi) * (phi*phi + 1) / (phi*phi - 1.)
                                           * (.5 * std::log(alpha * alphaMin / mappedAlphaZerothOrder / alphaMinMapped) - B);

              if (mappedNu > eps)
              {
                mappedAlphaFirstOrder = mappedNu * mappedNu / std::log(phi) * (phi*phi + 1) / (phi*phi - 1.)
                                           * (.5 * std::log(alpha * alphaMin / mappedAlphaZerothOrder / alphaMinMapped) - B);
              }
              else
              {
                mappedAlphaFirstOrder = nu * nu / std::log( arg ) / std::log( arg ) *
                                        (.5 * std::log(alpha * alphaMin / mappedAlphaZerothOrder / mappedAlphaZerothOrder ) - B);
              }

              return mappedAlphaZerothOrder * (1. + expiry * mappedAlphaFirstOrder);
            }
            else
            {
              double mappedAlphaZerothOrder = alpha;
              double mappedNu = effectiveZeroCorrelationNu(strike, forward, expiry);
              double mappedAlphaFirstOrder =
                          (1 - mappedNu*mappedNu/nu/nu - 1.5*rho*rho)*nu*nu/12. +
                          .25*beta*rho*alpha*nu*std::pow(forward, beta-1);
              return mappedAlphaZerothOrder * ( 1 + expiry * mappedAlphaFirstOrder );
            }
          }
        }
        double effectiveZeroCorrelationBeta(double strike,
                                            double forward,
                                            double expiry)  const override
        {
          return m_Beta->value(expiry);
        }
        double effectiveZeroCorrelationNu(double strike,
                                          double forward,
                                          double expiry)  const override
        {
          double alpha = m_Alpha->value(expiry);
          double beta = m_Beta->value(expiry);
          double rho = m_Rho->value(expiry);
          double nu = m_Nu->value(expiry);
          const double& eps = beagle::util::epsilon();

          if ( std::fabs( rho ) < eps )
            return nu;
          else
          {
            double nuSquared = nu * nu;
            double temp = nuSquared - 1.5*(nuSquared*rho*rho + alpha*nu*rho*(1-beta)*std::pow(forward, beta-1));
            if (temp < eps)
              return 0.;
            return std::sqrt( temp );
          }
        }
      protected:
        double kernel( double t, double s ) const
        {
          const double& pi = beagle::util::pi();
          if (!m_KernelMethod)
          {
            double factor = 2. * std::exp(-.125*t) / t / std::sqrt(pi*t);
            double coshS = std::cosh(s);

            real_func_t f = [=](double theta) {
              double tanTheta = std::tan(theta);
              double jacobian = 1 + tanTheta * tanTheta;
              double sPlusTanTheta = s + tanTheta;
              double exponentialTerm = std::exp(-.5 * sPlusTanTheta * sPlusTanTheta / t);
              double squareRootTerm = std::sqrt( std::cosh(sPlusTanTheta) - coshS );

              return jacobian * sPlusTanTheta * exponentialTerm * squareRootTerm;
            };
            return factor * m_QuadMethod->quadrature(beagle::math::RealFunction::createUnaryFunction(f), 0, .495 * pi );;
          }
          else
          {
            // When s is small enough, should replace the kernel with MacLaurin expansion to avoid infinity
            if ( std::fabs(s) < 5e-2 )
            {
              double factor = std::exp(-.5*s*s/t);
              double temp = std::exp(-.125*t);
              double expansion = 1 + ( 1./12 - (2688*t+80*t*t+21*t*t*t)/322560*temp)*s*s;
              return factor * expansion;
            }

            double sSquared = s*s;
            double sMinusSquared = 1. / sSquared;

            double g = s / std::tanh(s) - 1;
            double R = 1 + 3*t*g*.125*sMinusSquared
                     - 5*t*t*(-8*s*s+3*g*g+24*g)*sMinusSquared*sMinusSquared/128
                     + 35*t*t*t*(-40*sSquared+3*g*g*g+24*g*g+120*g)*sMinusSquared*sMinusSquared*sMinusSquared/1024;
            double deltaR = std::exp(.125*t) - (3072+384*t+24*t*t+t*t*t) / 3072;
            double factor = std::sqrt( std::sinh(s) / s ) * std::exp( -.5*sSquared/t - .125*t);
            return factor * (R + deltaR);
          }
        }
      private:
        virtual double doEffectiveStrike(double strike,
                                         double forward)  const = 0;
      protected:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        beagle::real_function_ptr_t m_Alpha;
        beagle::real_function_ptr_t m_Beta;
        beagle::real_function_ptr_t m_Rho;
        beagle::real_function_ptr_t m_Nu;
        bool m_KernelMethod;
        beagle::integration_method_ptr_t m_QuadMethod;
      };

      struct ClosedFormExactSABREuropeanOptionPricer : public ClosedFormMappedZeroCorrelationSABREuropeanOptionPricer
      {
        ClosedFormExactSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                const beagle::real_function_ptr_t& discounting,
                                                const beagle::real_function_ptr_t& alpha,
                                                const beagle::real_function_ptr_t& beta,
                                                const beagle::real_function_ptr_t& rho,
                                                const beagle::real_function_ptr_t& nu,
                                                bool useApproximateKernel,
                                                const beagle::integration_method_ptr_t& quadMethod) :
          ClosedFormMappedZeroCorrelationSABREuropeanOptionPricer(forward,
                                                                  discounting,
                                                                  alpha,
                                                                  beta,
                                                                  rho,
                                                                  nu,
                                                                  useApproximateKernel,
                                                                  quadMethod)
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
          const double& eps = beagle::util::epsilon();
          const double& pi = beagle::util::pi();

          double alphaToUse = effectiveZeroCorrelationAlpha(strike, forward, expiry);
          double betaToUse = effectiveZeroCorrelationBeta(strike, forward, expiry);
          double nuToUse = effectiveZeroCorrelationNu(strike, forward, expiry);

          if (nuToUse < eps)
          {
            return Pricer::formClosedFormExactCEVEuropeanOptionPricer(m_Forward,
                                                                      m_Discounting,
                                                                      beagle::math::RealFunction::createConstantFunction(betaToUse),
                                                                      beagle::math::RealFunction::createConstantFunction(alphaToUse),
                                                                      m_QuadMethod)->value(product);
          }

          double result = 0.;
          double oneMinusBeta = 1 - betaToUse;
          double eta = .5 / oneMinusBeta;
          double qK = std::pow( strike, oneMinusBeta ) / oneMinusBeta;
          double qF = std::pow( forward, oneMinusBeta ) / oneMinusBeta;
          double b = (qK * qK + qF * qF) / 2. / qK / qF;
          double tau = nuToUse * nuToUse * expiry;
          double vNaught = alphaToUse / nuToUse;
          double factor = std::sqrt( forward * strike ) / pi;

          real_func_t fOne = [=](double theta) {
            double temp = b - std::cos(theta);
            double s = std::asinh( std::sqrt(2*qK*qF*temp) / vNaught );

            return std::sin(eta * theta) * std::sin(theta) / temp
                   * kernel( tau, s ) / std::cosh(s);
          };
          result += m_QuadMethod->quadrature( beagle::math::RealFunction::createUnaryFunction(fOne), 0, pi );

          // double theta = 3e-4*PI;
          // double temp = b - std::cos(theta);
          // double s = std::asinh( std::sqrt(2*qK*qF*temp) / vNaught );
          // std::cout << s << "\t" << std::sin(eta * theta) * std::sin(theta) / temp << "\t" << kernel( tau, s ) << std::endl;

          real_func_t fTwo = [=](double theta) {
            double tanTheta = std::tan( theta );
            double expMinusX = std::exp( -tanTheta );
            double expMinusTwoX = expMinusX * expMinusX;
            double expMinusNuX = std::exp( -eta * tanTheta );
            double temp = 1 + expMinusTwoX + 2. * b * expMinusX;

            double tempOne = std::sqrt(( expMinusX + temp * qK * qF / vNaught / vNaught ) / expMinusX);
            return std::sin(eta * pi) * expMinusNuX * (1. - expMinusTwoX) / temp
                     * (1 + tanTheta * tanTheta)
                     * kernel( tau, std::acosh(tempOne) ) / tempOne;
          };
          result += m_QuadMethod->quadrature( beagle::math::RealFunction::createUnaryFunction(fTwo), 0, .495 * pi );

          result *= factor;
          if (pO->payoff()->isCall())
            result += std::max(forward - strike, 0.);
          else
            result += std::max(strike - forward, 0.);

          return result * discouting;
        }
        beagle::pricer_ptr_t updateModelParameters(const beagle::real_function_ptr_coll_t& params) const override
        {
          return Pricer::formClosedFormExactSABREuropeanOptionPricer(m_Forward,
                                                                     m_Discounting,
                                                                     params[0],
                                                                     params[1],
                                                                     params[2],
                                                                     params[3],
                                                                     m_KernelMethod,
                                                                     m_QuadMethod);
        }
      private:
        double doEffectiveStrike( double strike,
                                  double forward )  const override
        {
          return strike;
        }
      };

      struct ClosedFormFreeBoundarySABREuropeanOptionPricer : public ClosedFormMappedZeroCorrelationSABREuropeanOptionPricer
      {
        ClosedFormFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                       const beagle::real_function_ptr_t& discounting,
                                                       const beagle::real_function_ptr_t& alpha,
                                                       const beagle::real_function_ptr_t& beta,
                                                       const beagle::real_function_ptr_t& rho,
                                                       const beagle::real_function_ptr_t& nu,
                                                       bool useApproximateKernel,
                                                       const beagle::integration_method_ptr_t& quadMethod) :
          ClosedFormMappedZeroCorrelationSABREuropeanOptionPricer(forward,
                                                                  discounting,
                                                                  alpha,
                                                                  beta,
                                                                  rho,
                                                                  nu,
                                                                  useApproximateKernel,
                                                                  quadMethod)
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
          const double& eps = beagle::util::epsilon();
          const double& pi = beagle::util::pi();

          double alphaToUse = effectiveZeroCorrelationAlpha(strike, forward, expiry);
          double betaToUse = effectiveZeroCorrelationBeta(strike, forward, expiry);
          double nuToUse = effectiveZeroCorrelationNu(strike, forward, expiry);

          if (nuToUse < eps)
          {
            return Pricer::formClosedFormFreeBoundaryCEVEuropeanOptionPricer(m_Forward,
                                                                             m_Discounting,
                                                                             beagle::math::RealFunction::createConstantFunction(betaToUse),
                                                                             beagle::math::RealFunction::createConstantFunction(alphaToUse),
                                                                             m_QuadMethod)->value(product);
          }

          double result = 0.;
          double oneMinusBeta = 1 - betaToUse;
          double eta = .5 / oneMinusBeta;
          double qK = std::pow( std::fabs(strike), oneMinusBeta ) / oneMinusBeta;
          double qF = std::pow( std::fabs(forward), oneMinusBeta ) / oneMinusBeta;
          double b = (qK * qK + qF * qF) / 2. / qK / qF;
          double tau = nuToUse * nuToUse * expiry;
          double vNaught = nuToUse / alphaToUse;
          double factor = std::sqrt( std::fabs(forward * strike) ) / pi;

          real_func_t fOne = [=](double theta) {
            double temp = b - std::cos(theta);
            double s = std::asinh( std::sqrt(2*qK*qF*temp) * vNaught );

            return std::sin(eta * theta) * std::sin(theta) / temp
                   * kernel( tau, s ) / std::cosh(s);
          };
          if ( !(strike < eps) )
            result += m_QuadMethod->quadrature( beagle::math::RealFunction::createUnaryFunction(fOne), 0, pi );

          // double theta = 3e-4*PI;
          // double temp = b - std::cos(theta);
          // double s = std::asinh( std::sqrt(2*qK*qF*temp) / vNaught );
          // std::cout << s << "\t" << std::sin(eta * theta) * std::sin(theta) / temp << "\t" << kernel( tau, s ) << std::endl;

          real_func_t fTwo = [=](double theta) {
            double tanTheta = std::tan( theta );
            double expMinusX = std::exp( -tanTheta );
            double expMinusTwoX = expMinusX * expMinusX;
            double expMinusNuX = std::exp( -eta * tanTheta );
            double temp = 1 + expMinusTwoX + 2. * b * expMinusX;

            double tempOne = std::sqrt(( expMinusX + temp * qK * qF * vNaught * vNaught ) / expMinusX);
            double common = std::sin(eta * pi) * (1. - expMinusTwoX) / temp
                             * (1 + tanTheta * tanTheta)
                             * kernel( tau, std::acosh(tempOne) ) / tempOne;
            if ( !(strike < eps) )
              return common * std::cosh(eta * tanTheta);
            else
              return common * std::sinh(eta * tanTheta);
          };
          result += m_QuadMethod->quadrature( beagle::math::RealFunction::createUnaryFunction(fTwo), 0, .49 * pi );

          result *= factor;
          if (pO->payoff()->isCall())
            result += std::max(forward - strike, 0.);
          else
            result += std::max(strike - forward, 0.);

          return result * discouting;
        }
        beagle::pricer_ptr_t updateModelParameters( const beagle::real_function_ptr_coll_t& params ) const override
        {
          return Pricer::formClosedFormFreeBoundarySABREuropeanOptionPricer(m_Forward,
                                                                            m_Discounting,
                                                                            params[0],
                                                                            params[1],
                                                                            params[2],
                                                                            params[3],
                                                                            m_KernelMethod,
                                                                            m_QuadMethod);
        }
      private:
        double doEffectiveStrike( double strike,
                                  double forward )  const override
        {
          return std::max( strike, .1 * forward );
        }
      };
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormExactSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                        const beagle::real_function_ptr_t& discounting,
                                                        const beagle::real_function_ptr_t& alpha,
                                                        const beagle::real_function_ptr_t& beta,
                                                        const beagle::real_function_ptr_t& rho,
                                                        const beagle::real_function_ptr_t& nu,
                                                        bool useApproximateKernel,
                                                        const beagle::integration_method_ptr_t& quadMethod)
    {
      return std::make_shared<impl::ClosedFormExactSABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, useApproximateKernel, quadMethod);
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                               const beagle::real_function_ptr_t& discounting,
                                                               const beagle::real_function_ptr_t& alpha,
                                                               const beagle::real_function_ptr_t& beta,
                                                               const beagle::real_function_ptr_t& rho,
                                                               const beagle::real_function_ptr_t& nu,
                                                               bool useApproximateKernel,
                                                               const beagle::integration_method_ptr_t& quadMethod)
    {
      return std::make_shared<impl::ClosedFormFreeBoundarySABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, useApproximateKernel, quadMethod);
    }
  }
}
