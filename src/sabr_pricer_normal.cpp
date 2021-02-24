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
      struct ClsoedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer : public ClosedFormEuropeanOptionPricer
      {
        ClsoedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                const beagle::real_function_ptr_t& discounting,
                                                                const beagle::real_function_ptr_t& alpha,
                                                                const beagle::real_function_ptr_t& beta,
                                                                const beagle::real_function_ptr_t& rho,
                                                                const beagle::real_function_ptr_t& nu,
                                                                bool useApproximateKernel,
                                                                const beagle::integration_method_ptr_t& quadMethod) :
          ClosedFormEuropeanOptionPricer(forward, discounting),
          m_Alpha(alpha),
          m_Beta(beta),
          m_Rho(rho),
          m_Nu(nu),
          m_KernelMethod(useApproximateKernel),
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
          double alpha = m_Alpha->value(expiry);
          double beta = m_Beta->value(expiry);
          double forward = forwardCurve()->value(expiry);
          double effectiveAlpha = alpha * std::pow(forward, beta);

          beagle::pricer_ptr_t pFreeBoundarySABRPricer = Pricer::formClosedFormFreeBoundarySABREuropeanOptionPricer(forwardCurve(),
                                                                                                                    discountCurve(),
                                                                                                                    m_Alpha,
                                                                                                                    m_Beta,
                                                                                                                    m_Rho,
                                                                                                                    m_Nu,
                                                                                                                    m_KernelMethod,
                                                                                                                    m_QuadMethod);
          beagle::pricer_ptr_t pNormalFreeBoundarySABRPricer = Pricer::formClosedFormNormalFreeBoundarySABREuropeanOptionPricer(forwardCurve(),
                                                                                                                                discountCurve(),
                                                                                                                                beagle::math::RealFunction::createConstantFunction(effectiveAlpha),
                                                                                                                                m_Rho,
                                                                                                                                m_Nu,
                                                                                                                                m_KernelMethod,
                                                                                                                                m_QuadMethod);

          beagle::pricer_ptr_t pControlFreeBoundarySABRPricer = Pricer::formClosedFormFreeBoundarySABREuropeanOptionPricer(forwardCurve(),
                                                                                                                           discountCurve(),
                                                                                                                           beagle::math::RealFunction::createConstantFunction(effectiveAlpha),
                                                                                                                           beagle::math::RealFunction::createConstantFunction(beagle::util::epsilon()),
                                                                                                                           m_Rho,
                                                                                                                           m_Nu,
                                                                                                                           m_KernelMethod,
                                                                                                                           m_QuadMethod);

          return pFreeBoundarySABRPricer->value(product)
                  + pNormalFreeBoundarySABRPricer->value(product)
                  - pControlFreeBoundarySABRPricer->value(product);
        }
        beagle::real_function_ptr_coll_t modelParameters( void ) const override
        {
          return {m_Alpha, m_Beta, m_Rho, m_Nu};
        }
        beagle::pricer_ptr_t updateModelParameters( const beagle::real_function_ptr_coll_t& params ) const override
        {
          return Pricer::formClosedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer(forwardCurve(),
                                                                                          discountCurve(),
                                                                                          params[0],
                                                                                          params[1],
                                                                                          params[2],
                                                                                          params[3],
                                                                                          m_KernelMethod,
                                                                                          m_QuadMethod);
        }
        int numberOfParameters( void ) const override
        {
          return 4;
        }
      private:
        beagle::real_function_ptr_t m_Alpha;
        beagle::real_function_ptr_t m_Beta;
        beagle::real_function_ptr_t m_Rho;
        beagle::real_function_ptr_t m_Nu;
        bool m_KernelMethod;
        beagle::integration_method_ptr_t m_QuadMethod;
      };

      struct ClsoedFormNormalFreeBoundarySABREuropeanOptionPricer : public ClosedFormEuropeanOptionPricer
      {
        ClsoedFormNormalFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                             const beagle::real_function_ptr_t& discounting,
                                                             const beagle::real_function_ptr_t& alpha,
                                                             const beagle::real_function_ptr_t& rho,
                                                             const beagle::real_function_ptr_t& nu,
                                                             bool useApproximateKernel,
                                                             const beagle::integration_method_ptr_t& quadMethod) :
          ClosedFormEuropeanOptionPricer(forward, discounting),
          m_Alpha(alpha),
          m_Rho(rho),
          m_Nu(nu),
          m_KernelMethod(useApproximateKernel),
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
          double forward = forwardCurve()->value(expiry);
          double discounting = discountCurve()->value(expiry);
          const double& pi = beagle::util::pi();

          double alpha = m_Alpha->value(expiry);
          double rho = m_Rho->value(expiry);
          double nu = m_Nu->value(expiry);
          double vNaught = alpha / nu;
          double k = (strike - forward) / vNaught + rho;
          double sNaught = std::acosh( (std::sqrt(k*k+1.-rho*rho) - rho*k) / (1. - rho*rho) );
          double factor = vNaught / pi;

          real_func_t fOne = [=](double theta) {
            double s = std::tan(theta);
            double jacobian = 1. + s*s;
            double ker = kernel(nu*nu*expiry, sNaught + s);
            double denominator = std::sinh(sNaught + s);
            double temp = (k - rho*std::cosh(sNaught + s)) / denominator;

            return ker * jacobian * std::sqrt(1. - temp*temp);
          };

          double result = factor * m_QuadMethod->quadrature( beagle::math::RealFunction::createUnaryFunction(fOne), 0, .49*pi );

          if (pO->payoff()->isCall())
            result += std::max(forward - strike, 0.);
          else
            result += std::max(strike - forward, 0.);

          return result * discounting;
        }
        beagle::real_function_ptr_coll_t modelParameters( void ) const override
        {
          return {m_Alpha, m_Rho, m_Nu};
        }
        beagle::pricer_ptr_t updateModelParameters( const beagle::real_function_ptr_coll_t& params ) const override
        {
          return Pricer::formClosedFormNormalFreeBoundarySABREuropeanOptionPricer(forwardCurve(),
                                                                                  discountCurve(),
                                                                                  params[0],
                                                                                  params[1],
                                                                                  params[2],
                                                                                  m_KernelMethod,
                                                                                  m_QuadMethod);
        }
        int numberOfParameters( void ) const override
        {
          return 3;
        }
      private:
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
            return factor * m_QuadMethod->quadrature( beagle::math::RealFunction::createUnaryFunction(f), 0, .49 * pi );;
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
        beagle::real_function_ptr_t m_Alpha;
        beagle::real_function_ptr_t m_Rho;
        beagle::real_function_ptr_t m_Nu;
        bool m_KernelMethod;
        beagle::integration_method_ptr_t m_QuadMethod;
      };
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                             const beagle::real_function_ptr_t& discounting,
                                                                             const beagle::real_function_ptr_t& alpha,
                                                                             const beagle::real_function_ptr_t& beta,
                                                                             const beagle::real_function_ptr_t& rho,
                                                                             const beagle::real_function_ptr_t& nu,
                                                                             bool useApproximateKernel,
                                                                             const beagle::integration_method_ptr_t& quadMethod)
    {
      return std::make_shared<impl::ClsoedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, useApproximateKernel, quadMethod);
    }

    beagle::pricer_ptr_t
    Pricer::formClosedFormNormalFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                                                     const beagle::real_function_ptr_t& discounting,
                                                                     const beagle::real_function_ptr_t& alpha,
                                                                     const beagle::real_function_ptr_t& rho,
                                                                     const beagle::real_function_ptr_t& nu,
                                                                     bool useApproximateKernel,
                                                                     const beagle::integration_method_ptr_t& quadMethod)
    {
      return std::make_shared<impl::ClsoedFormNormalFreeBoundarySABREuropeanOptionPricer>(forward, discounting, alpha, rho, nu, useApproximateKernel, quadMethod);
    }
  }
}
