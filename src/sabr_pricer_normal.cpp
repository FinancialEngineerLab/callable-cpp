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
      // struct NormalImprovedFreeBoundarySABR : public OptionPricer
      // {
      //   NormalImprovedFreeBoundarySABR( double alpha,
      //                                   double beta,
      //                                   double rho,
      //                                   double nu,
      //                                   bool useApproximateKernel,
      //                                   const integration_ptr_t& quadMethod ) :
      //     m_Alpha( alpha ),
      //     m_Beta( beta ),
      //     m_Rho( rho ),
      //     m_Nu( nu ),
      //     m_Kernel( useApproximateKernel ),
      //     m_QuadMethod( quadMethod )
      //   {
      //     if (alpha < EPSILON)
      //     {
      //       std::cout << "alpha = " << alpha << std::endl;
      //       abort();
      //     }
      //     if (beta < EPSILON || beta-.5 > EPSILON)
      //     {
      //       std::cout << "beta = " << beta << std::endl;
      //       abort();
      //     }
      //     if (rho+1. < EPSILON || rho-1. > EPSILON)
      //     {
      //       std::cout << "rho = " << rho << std::endl;
      //       abort();
      //     }
      //     if (nu < EPSILON)
      //     {
      //       std::cout << "nu = " << nu << std::endl;
      //       abort();
      //     }
      //   }
      //   virtual ~NormalImprovedFreeBoundarySABR( void ) { }
      // public:
      //   const std::string& description( void ) const override
      //   {
      //     static std::string name("NormalImprovedFreeBoundarySABR");
      //     return name;
      //   }
      //   double callValue( double strike,
      //                     double forward,
      //                     double expiry ) const override
      //   {
      //     opt_ptr_t pFreeBoundarySABRPricer = OptionPricer::freeBoundarySABR(m_Alpha, m_Beta, m_Rho, m_Nu, m_Kernel, m_QuadMethod);

      //     double effectiveAlpha = m_Alpha * std::pow(forward, m_Beta);
      //     opt_ptr_t pNormalFreeBoundarySABRPricer = OptionPricer::normalFreeBoundarySABR(effectiveAlpha, m_Rho, m_Nu, m_Kernel, m_QuadMethod);
      //     opt_ptr_t pControlFreeBoundarySABRPricer = OptionPricer::freeBoundarySABR(effectiveAlpha, EPSILON, m_Rho, m_Nu, m_Kernel, m_QuadMethod);

      //     return std::max(forward - strike, 0.)
      //             + pFreeBoundarySABRPricer->callValue( strike, forward, expiry )
      //             + pNormalFreeBoundarySABRPricer->callValue( strike, forward, expiry )
      //             - pControlFreeBoundarySABRPricer->callValue( strike, forward, expiry );
      //   }
      //   dbl_vec_t modelParameters( void ) const override
      //   {
      //     dbl_vec_t params(4);
      //     params << m_Alpha, m_Beta, m_Rho, m_Nu;

      //     return params;
      //   }
      //   opt_ptr_t updateModelParameters( const dbl_vec_t& params ) const override
      //   {
      //     return std::make_shared<NormalImprovedFreeBoundarySABR>(params(0), params(1), params(2), params(3), m_Kernel, m_QuadMethod);
      //   }
      //   int numberOfRiskFactors( void ) const override
      //   {
      //     return 4;
      //   }
      // private:
      //   double m_Alpha;
      //   double m_Beta;
      //   double m_Rho;
      //   double m_Nu;
      //   bool m_Kernel;
      //   integration_ptr_t m_QuadMethod;
      // };

      // struct NormalFreeBoundarySABR : public OptionPricer
      // {
      //   NormalFreeBoundarySABR( double alpha,
      //                           double rho,
      //                           double nu,
      //                           bool useApproximateKernel,
      //                           const integration_ptr_t& quadMethod ) :
      //     m_Alpha( alpha ),
      //     m_Rho( rho ),
      //     m_Nu( nu ),
      //     m_Kernel( useApproximateKernel ),
      //     m_QuadMethod( quadMethod )
      //   {
      //     if (alpha < EPSILON)
      //     {
      //       std::cout << "alpha = " << alpha << std::endl;
      //       abort();
      //     }
      //     if (rho+1. < EPSILON || rho-1. > EPSILON)
      //     {
      //       std::cout << "rho = " << rho << std::endl;
      //       abort();
      //     }
      //     if (nu < EPSILON)
      //     {
      //       std::cout << "nu = " << nu << std::endl;
      //       abort();
      //     }
      //   }
      //   virtual ~NormalFreeBoundarySABR( void ) { }
      // public:
      //   const std::string& description( void ) const override
      //   {
      //     static std::string name(std::string("NormalFreeBoundarySABR:Kernel") + (m_Kernel ? "Series:" : "Integral:") + std::to_string(m_QuadMethod->numberOfSteps()));
      //     return name;
      //   }
      //   double callValue( double strike,
      //                     double forward,
      //                     double expiry ) const override
      //   {
      //     double vNaught = m_Alpha / m_Nu;
      //     double k = (strike - forward) / vNaught + m_Rho;
      //     double sNaught = std::acosh( (std::sqrt(k*k+1.-m_Rho*m_Rho) - m_Rho*k) / (1. - m_Rho*m_Rho) );
      //     double factor = vNaught / PI;

      //     real_func_t fOne = [=](double theta) {
      //       double s = std::tan(theta);
      //       double jacobian = 1. + s*s;
      //       double ker = kernel(m_Nu*m_Nu*expiry, sNaught + s);
      //       double denominator = std::sinh(sNaught + s);
      //       double temp = (k - m_Rho*std::cosh(sNaught + s)) / denominator;

      //       return ker * jacobian * std::sqrt(1. - temp*temp);
      //     };

      //     return factor * m_QuadMethod->quadrature( fOne, 0, .49*PI ); // + std::max(forward - strike, 0.);
      //   }
      //   dbl_vec_t modelParameters( void ) const override
      //   {
      //     dbl_vec_t params(3);
      //     params << m_Alpha, m_Rho, m_Nu;

      //     return params;
      //   }
      //   opt_ptr_t updateModelParameters( const dbl_vec_t& params ) const override
      //   {
      //     return std::make_shared<NormalFreeBoundarySABR>(params(0), params(1), params(2), m_Kernel, m_QuadMethod);
      //   }
      //   int numberOfRiskFactors( void ) const override
      //   {
      //     return 3;
      //   }
      // private:
      //   double kernel( double t, double s ) const
      //   {
      //     if (!m_Kernel)
      //     {
      //       double factor = 2. * std::exp(-.125*t) / t / std::sqrt(PI*t);
      //       double coshS = std::cosh(s);

      //       real_func_t f = [=](double theta) {
      //         double tanTheta = std::tan(theta);
      //         double jacobian = 1 + tanTheta * tanTheta;
      //         double sPlusTanTheta = s + tanTheta;
      //         double exponentialTerm = std::exp(-.5 * sPlusTanTheta * sPlusTanTheta / t);
      //         double squareRootTerm = std::sqrt( std::cosh(sPlusTanTheta) - coshS );

      //         return jacobian * sPlusTanTheta * exponentialTerm * squareRootTerm;
      //       };
      //       return factor * m_QuadMethod->quadrature( f, 0, .49 * PI );;
      //     }
      //     else
      //     {
      //       // When s is small enough, should replace the kernel with MacLaurin expansion to avoid infinity
      //       if ( std::fabs(s) < 5e-2 )
      //       {
      //         double factor = std::exp(-.5*s*s/t);
      //         double temp = std::exp(-.125*t);
      //         double expansion = 1 + ( 1./12 - (2688*t+80*t*t+21*t*t*t)/322560*temp)*s*s;
      //         return factor * expansion;
      //       }

      //       double sSquared = s*s;
      //       double sMinusSquared = 1. / sSquared;

      //       double g = s / std::tanh(s) - 1;
      //       double R = 1 + 3*t*g*.125*sMinusSquared
      //                - 5*t*t*(-8*s*s+3*g*g+24*g)*sMinusSquared*sMinusSquared/128
      //                + 35*t*t*t*(-40*sSquared+3*g*g*g+24*g*g+120*g)*sMinusSquared*sMinusSquared*sMinusSquared/1024;
      //       double deltaR = std::exp(.125*t) - (3072+384*t+24*t*t+t*t*t) / 3072;
      //       double factor = std::sqrt( std::sinh(s) / s ) * std::exp( -.5*sSquared/t - .125*t);
      //       return factor * (R + deltaR);
      //     }
      //   }
      // private:
      //   double m_Alpha;
      //   double m_Rho;
      //   double m_Nu;
      //   bool m_Kernel;
      //   integration_ptr_t m_QuadMethod;
      // };
    }

    // beagle::pricer_ptr_t formClosedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
    //                                                                                       const beagle::real_function_ptr_t& discounting,
    //                                                                                       const beagle::real_function_ptr_t& alpha,
    //                                                                                       const beagle::real_function_ptr_t& beta,
    //                                                                                       const beagle::real_function_ptr_t& rho,
    //                                                                                       const beagle::real_function_ptr_t& nu,
    //                                                                                       bool useApproximateKernel,
    //                                                                                       const beagle::integration_method_ptr_t& quadMethod)
    // {
    //   return std::make_shared<impl::ClsoedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, useApproximateKernel, quadMethod);
    // }

    // beagle::pricer_ptr_t formClosedFormNormalFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
    //                                                                               const beagle::real_function_ptr_t& discounting,
    //                                                                               const beagle::real_function_ptr_t& alpha,
    //                                                                               const beagle::real_function_ptr_t& beta,
    //                                                                               const beagle::real_function_ptr_t& rho,
    //                                                                               const beagle::real_function_ptr_t& nu,
    //                                                                               bool useApproximateKernel,
    //                                                                               const beagle::integration_method_ptr_t& quadMethod)
    // {
    //   return std::make_shared<impl::ClsoedFormNormalImprovedFreeBoundarySABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, useApproximateKernel, quadMethod);
    // }
  }
}
