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
      namespace
      {
        void checkSABRParameters(double alpha,
                                 double beta,
                                 double rho,
                                 double nu,
                                 bool isFreeBoundarySABR)
        {
          double betaMax = isFreeBoundarySABR ? .5 : 1.;
          const double& eps = util::epsilon();
          if (alpha < eps)
            throw("The alpha parameter of the SABR model must be positive!");
          if (beta < eps || beta - betaMax > eps)
            throw("The beta parameter of the SABR model must be between zero and " + std::to_string(betaMax) + "!");
          if (rho + 1. < eps || rho - 1. > eps)
            throw("The rho parameter of the SABR model must be between minus one and one!");
          if (nu < eps)
            throw("The nu parameter of the SABR model must be positive!");
        }
      }

      struct ClosedFormSABREuropeanOptionPricer : public ClosedFormEuropeanOptionPricer
      {
        ClosedFormSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
                                           const beagle::real_function_ptr_t& discounting,
                                           const beagle::real_function_ptr_t& alpha,
                                           const beagle::real_function_ptr_t& beta,
                                           const beagle::real_function_ptr_t& rho,
                                           const beagle::real_function_ptr_t& nu) :
          m_Forward(forward),
          m_Discounting(discounting),
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
          double forward = m_Forward->value(expiry);
          double discouting = m_Discounting->value(expiry);

          double result = util::bsCall(strike, forward, expiry, impliedBlackScholesVolatility(product));
          const auto& payoff = pO->payoff();
          if (payoff->isPut())
            result -= (forward - strike);

          return result * discouting;
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

          double forward = m_Forward->value(expiry);
          double alpha = m_Alpha->value(expiry);
          double beta = m_Beta->value(expiry);
          double rho = m_Rho->value(expiry);
          double nu = m_Nu->value(expiry);
          checkSABRParameters(alpha, beta, rho, nu, false);

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

          if ( std::fabs( forward - strike ) > util::epsilon() )
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
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
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
          m_Forward(forward),
          m_Discounting(discounting),
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
          double forward = m_Forward->value(expiry);
          double discouting = m_Discounting->value(expiry);

          double result = util::bsCall(strike + m_Shift, forward + m_Shift, expiry, impliedBlackScholesVolatility(product));
          const auto& payoff = pO->payoff();
          if (payoff->isPut())
            result -= (forward - strike);

          return result * discouting;
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

          double forward = m_Forward->value(expiry) + m_Shift;
          double alpha = m_Alpha->value(expiry);
          double beta = m_Beta->value(expiry);
          double rho = m_Rho->value(expiry);
          double nu = m_Nu->value(expiry);
          checkSABRParameters(alpha, beta, rho, nu, false);

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

          if ( std::fabs( forward - strike ) > util::epsilon() )
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
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        beagle::real_function_ptr_t m_Alpha;
        beagle::real_function_ptr_t m_Beta;
        beagle::real_function_ptr_t m_Rho;
        beagle::real_function_ptr_t m_Nu;
        double m_Shift;
      };

      // struct ZeroCorrelationMappedSABR : public OptionPricer,
      //                                    public mixins::MappedZeroCorrelationSABRParameters
      // {
      //   ZeroCorrelationMappedSABR( double alpha,
      //                              double beta,
      //                              double rho,
      //                              double nu,
      //                              bool useApproximateKernel,
      //                              const integration_ptr_t& quadMethod ) :
      //     m_Alpha( alpha ),
      //     m_Beta( beta ),
      //     m_Rho( rho ),
      //     m_Nu( nu ),
      //     m_Kernel( useApproximateKernel ),
      //     m_QuadMethod( quadMethod )
      //   {
      //   }
      //   virtual ~ZeroCorrelationMappedSABR( void ) { }
      // public:
      //   virtual const std::string& description( void ) const = 0;
      //   virtual double callValue( double strike,
      //                             double forward,
      //                             double expiry ) const = 0;
      //   virtual opt_ptr_t updateModelParameters( const dbl_vec_t& params ) const = 0;

      //   dbl_vec_t modelParameters( void ) const override
      //   {
      //     dbl_vec_t params(4);
      //     params << m_Alpha, m_Beta, m_Rho, m_Nu;

      //     return params;
      //   }
      //   int numberOfRiskFactors( void ) const override
      //   {
      //     return 4;
      //   }
      //   double effectiveZeroCorrelationAlpha( double strike, double forward, double expiry )  const override
      //   {
      //     if ( std::fabs( m_Rho ) < EPSILON )
      //       return m_Alpha;
      //     else
      //     {
      //       double mappedAlphaZerothOrder;
      //       double mappedAlphaFirstOrder;
      //       if ( std::fabs( forward - strike ) > EPSILON )
      //       {
      //         double mappedBeta = effectiveZeroCorrelationBeta();
      //         double mappedNu = effectiveZeroCorrelationNu( forward );

      //         double effectiveStrike = std::fabs(doEffectiveStrike( strike, forward ));
      //         double effectiveForward = std::fabs( forward );

      //         double oneMinusBeta = 1 - m_Beta;
      //         double qK = std::pow( effectiveStrike, oneMinusBeta ) / oneMinusBeta;
      //         double qF = std::pow( effectiveForward, oneMinusBeta ) / oneMinusBeta;
      //         double deltaQ = qK - qF;

      //         double oneMinusMappedBeta = 1. - mappedBeta;
      //         double qKMapped = std::pow( effectiveStrike, oneMinusMappedBeta ) / oneMinusMappedBeta;
      //         double qFMapped = std::pow( effectiveForward, oneMinusMappedBeta ) / oneMinusMappedBeta;
      //         double deltaQMapped = qKMapped - qFMapped;

      //         double alphaMinSquared = m_Nu*m_Nu*deltaQ*deltaQ + 2*m_Rho*m_Nu*deltaQ*m_Alpha + m_Alpha*m_Alpha;
      //         double alphaMin = std::sqrt( alphaMinSquared );
      //         double arg = (alphaMin + m_Rho*m_Alpha + m_Nu*deltaQ) / (1.+m_Rho) / m_Alpha;
      //         double phi = std::pow( arg, mappedNu/m_Nu );

      //         if (mappedNu < EPSILON)
      //           mappedAlphaZerothOrder = m_Nu * deltaQMapped / std::log( arg );
      //         else
      //           mappedAlphaZerothOrder = 2.*phi*deltaQMapped*mappedNu / (phi*phi - 1);

      //         double alphaMinMapped =
      //               std::sqrt(mappedNu*mappedNu*deltaQMapped*deltaQMapped + mappedAlphaZerothOrder*mappedAlphaZerothOrder);
      //         double rootOneMinusRhoSquared = std::sqrt(1. - m_Rho*m_Rho);
      //         double L = alphaMin / qK / m_Nu / rootOneMinusRhoSquared;
      //         double uNaught = (deltaQ * m_Nu * m_Rho + m_Alpha - alphaMin) / deltaQ / m_Nu / rootOneMinusRhoSquared;

      //         double I = 2. * uNaught / (1 + uNaught);
      //         if (L - 1. < 0)
      //         {
      //           double rootOneMinusLSquared = std::sqrt(1. - L * L);
      //           double tempOne = (uNaught + L) / rootOneMinusLSquared;
      //           double tempTwo = L / rootOneMinusLSquared;
      //           I = 2. * (std::atan(tempOne) - std::atan(tempTwo)) / rootOneMinusLSquared;
      //         }
      //         else if (L - 1. > 0)
      //         {
      //           double rootLSquaredMinusOne = std::sqrt(L * L - 1.);
      //           double tempOne = L + rootLSquaredMinusOne;
      //           double tempTwo = L - rootLSquaredMinusOne;
      //           I = (std::log(1. + uNaught * tempOne) - std::log(1. + uNaught * tempTwo)) / rootLSquaredMinusOne;
      //         }

      //         double B = -.5 * m_Beta / oneMinusBeta * m_Rho / rootOneMinusRhoSquared
      //                  * (PI - std::acos( -(deltaQ * m_Nu + m_Alpha * m_Rho) / alphaMin ) - std::acos( m_Rho ) - I);

      //         double mappedAlphaFirstOrder = mappedNu * mappedNu / std::log(phi) * (phi*phi + 1) / (phi*phi - 1.)
      //                                      * (.5 * std::log(m_Alpha * alphaMin / mappedAlphaZerothOrder / alphaMinMapped) - B);

      //         if (mappedNu > EPSILON)
      //         {
      //           mappedAlphaFirstOrder = mappedNu * mappedNu / std::log(phi) * (phi*phi + 1) / (phi*phi - 1.)
      //                                      * (.5 * std::log(m_Alpha * alphaMin / mappedAlphaZerothOrder / alphaMinMapped) - B);
      //         }
      //         else
      //         {
      //           mappedAlphaFirstOrder = m_Nu * m_Nu / std::log( arg ) / std::log( arg ) *
      //                                   (.5 * std::log(m_Alpha * alphaMin / mappedAlphaZerothOrder / mappedAlphaZerothOrder ) - B);
      //         }

      //         return mappedAlphaZerothOrder * (1. + expiry * mappedAlphaFirstOrder);
      //       }
      //       else
      //       {
      //         double mappedAlphaZerothOrder = m_Alpha;
      //         double mappedNu = effectiveZeroCorrelationNu( forward );
      //         double mappedAlphaFirstOrder =
      //                     (1 - mappedNu*mappedNu/m_Nu/m_Nu - 1.5*m_Rho*m_Rho)*m_Nu*m_Nu/12. +
      //                     .25*m_Beta*m_Rho*m_Alpha*m_Nu*std::pow(forward, m_Beta-1);
      //         return mappedAlphaZerothOrder * ( 1 + expiry * mappedAlphaFirstOrder );
      //       }
      //     }
      //   }
      //   double effectiveZeroCorrelationBeta( void )  const override
      //   {
      //     return m_Beta;
      //   }
      //   double effectiveZeroCorrelationNu( double forward )  const override
      //   {
      //     if ( std::fabs( m_Rho ) < EPSILON )
      //       return m_Nu;
      //     else
      //     {
      //       double nuSquared = m_Nu * m_Nu;
      //       double temp = nuSquared - 1.5*(nuSquared*m_Rho*m_Rho + m_Alpha*m_Nu*m_Rho*(1-m_Beta)*std::pow(forward, m_Beta-1));
      //       if (temp < EPSILON)
      //         return 0.;
      //       return std::sqrt( temp );
      //     }
      //   }
      // protected:
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
      //   bool kernelMethod( void ) const
      //   {
      //     return m_Kernel;
      //   }
      //   const integration_ptr_t& integrationMethod( void ) const
      //   {
      //     return m_QuadMethod;
      //   }
      // private:
      //   virtual double doEffectiveStrike( double strike,
      //                                     double forward )  const = 0;
      // private:
      //   double m_Alpha;
      //   double m_Beta;
      //   double m_Rho;
      //   double m_Nu;
      //   bool m_Kernel;
      //   integration_ptr_t m_QuadMethod;
      // };

      // struct ExactSABR : public ZeroCorrelationMappedSABR
      // {
      //   ExactSABR( double alpha,
      //              double beta,
      //              double rho,
      //              double nu,
      //              bool useApproximateKernel,
      //              const integration_ptr_t& quadMethod ) :
      //     ZeroCorrelationMappedSABR( alpha, beta, rho, nu, useApproximateKernel, quadMethod )
      //   {
      //     if (alpha < EPSILON)
      //     {
      //       std::cout << "alpha = " << alpha << std::endl;
      //       abort();
      //     }
      //     if (beta < EPSILON || beta-1. > EPSILON)
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
      //   virtual ~ExactSABR( void ) { }
      // public:
      //   const std::string& description( void ) const override
      //   {
      //     static std::string name(std::string("ExactSABR:Kernel") + (kernelMethod() ? "Series:" : "Integral:") + std::to_string(integrationMethod()->numberOfSteps()));
      //     return name;
      //   }
      //   double callValue( double strike,
      //                     double forward,
      //                     double expiry ) const override
      //   {
      //     double alphaToUse = effectiveZeroCorrelationAlpha( strike, forward, expiry );
      //     double betaToUse = effectiveZeroCorrelationBeta();
      //     double nuToUse = effectiveZeroCorrelationNu( forward );

      //     if (nuToUse < EPSILON)
      //     {
      //       return OptionPricer::exactCEV(betaToUse, alphaToUse, integrationMethod())->callValue(strike, forward, expiry);
      //     }

      //     double result = 0.;
      //     double oneMinusBeta = 1 - betaToUse;
      //     double eta = .5 / oneMinusBeta;
      //     double qK = std::pow( strike, oneMinusBeta ) / oneMinusBeta;
      //     double qF = std::pow( forward, oneMinusBeta ) / oneMinusBeta;
      //     double b = (qK * qK + qF * qF) / 2. / qK / qF;
      //     double tau = nuToUse * nuToUse * expiry;
      //     double vNaught = alphaToUse / nuToUse;
      //     double factor = std::sqrt( forward * strike ) / PI;

      //     real_func_t fOne = [=](double theta) {
      //       double temp = b - std::cos(theta);
      //       double s = std::asinh( std::sqrt(2*qK*qF*temp) / vNaught );

      //       return std::sin(eta * theta) * std::sin(theta) / temp
      //              * kernel( tau, s ) / std::cosh(s);
      //     };
      //     result += integrationMethod()->quadrature( fOne, 0, PI );

      //     // double theta = 3e-4*PI;
      //     // double temp = b - std::cos(theta);
      //     // double s = std::asinh( std::sqrt(2*qK*qF*temp) / vNaught );
      //     // std::cout << s << "\t" << std::sin(eta * theta) * std::sin(theta) / temp << "\t" << kernel( tau, s ) << std::endl;

      //     real_func_t fTwo = [=](double theta) {
      //       double tanTheta = std::tan( theta );
      //       double expMinusX = std::exp( -tanTheta );
      //       double expMinusTwoX = expMinusX * expMinusX;
      //       double expMinusNuX = std::exp( -eta * tanTheta );
      //       double temp = 1 + expMinusTwoX + 2. * b * expMinusX;

      //       double tempOne = std::sqrt(( expMinusX + temp * qK * qF / vNaught / vNaught ) / expMinusX);
      //       return std::sin(eta * PI) * expMinusNuX * (1. - expMinusTwoX) / temp
      //                * (1 + tanTheta * tanTheta)
      //                * kernel( tau, std::acosh(tempOne) ) / tempOne;
      //     };
      //     result += integrationMethod()->quadrature( fTwo, 0, .495 * PI );
      //     return factor * result + std::max(forward - strike, 0.);
      //   }
      //   opt_ptr_t updateModelParameters( const dbl_vec_t& params ) const override
      //   {
      //     return std::make_shared<ExactSABR>(params(0), params(1), params(2), params(3), kernelMethod(), integrationMethod());
      //   }
      // private:
      //   double doEffectiveStrike( double strike,
      //                             double forward )  const override
      //   {
      //     return strike;
      //   }
      // };

      // struct FreeBoundarySABR : public ZeroCorrelationMappedSABR
      // {
      //   FreeBoundarySABR( double alpha,
      //                     double beta,
      //                     double rho,
      //                     double nu,
      //                     bool useApproximateKernel,
      //                     const integration_ptr_t& quadMethod ) :
      //     ZeroCorrelationMappedSABR( alpha, beta, rho, nu, useApproximateKernel, quadMethod )
      //   {
      //     if (alpha < EPSILON)
      //     {
      //       std::cout << "alpha = " << alpha << std::endl;
      //       abort();
      //     }
      //     if (beta < 0. || beta-.5 > 0.)
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
      //   virtual ~FreeBoundarySABR( void ) { }
      // public:
      //   const std::string& description( void ) const override
      //   {
      //     static std::string name(std::string("FreeBoundarySABR:Kernel") + (kernelMethod() ? "Series:" : "Integral:") + std::to_string(integrationMethod()->numberOfSteps()));
      //     return name;
      //   }
      //   double callValue( double strike,
      //                     double forward,
      //                     double expiry ) const override
      //   {
      //     double alphaToUse = effectiveZeroCorrelationAlpha( strike, forward, expiry );
      //     double betaToUse = effectiveZeroCorrelationBeta();
      //     double nuToUse = effectiveZeroCorrelationNu( forward );

      //     if (nuToUse < EPSILON)
      //     {
      //       return OptionPricer::freeBoundaryCEV(betaToUse, alphaToUse, integrationMethod())->callValue(strike, forward, expiry);
      //     }

      //     double result = 0.;
      //     double oneMinusBeta = 1 - betaToUse;
      //     double eta = .5 / oneMinusBeta;
      //     double qK = std::pow( std::fabs(strike), oneMinusBeta ) / oneMinusBeta;
      //     double qF = std::pow( std::fabs(forward), oneMinusBeta ) / oneMinusBeta;
      //     double b = (qK * qK + qF * qF) / 2. / qK / qF;
      //     double tau = nuToUse * nuToUse * expiry;
      //     double vNaught = nuToUse / alphaToUse;
      //     double factor = std::sqrt( std::fabs(forward * strike) ) / PI;

      //     real_func_t fOne = [=](double theta) {
      //       double temp = b - std::cos(theta);
      //       double s = std::asinh( std::sqrt(2*qK*qF*temp) * vNaught );

      //       return std::sin(eta * theta) * std::sin(theta) / temp
      //              * kernel( tau, s ) / std::cosh(s);
      //     };
      //     if ( !(strike < EPSILON) )
      //       result += integrationMethod()->quadrature( fOne, 0, PI );

      //     // double theta = 3e-4*PI;
      //     // double temp = b - std::cos(theta);
      //     // double s = std::asinh( std::sqrt(2*qK*qF*temp) / vNaught );
      //     // std::cout << s << "\t" << std::sin(eta * theta) * std::sin(theta) / temp << "\t" << kernel( tau, s ) << std::endl;

      //     real_func_t fTwo = [=](double theta) {
      //       double tanTheta = std::tan( theta );
      //       double expMinusX = std::exp( -tanTheta );
      //       double expMinusTwoX = expMinusX * expMinusX;
      //       double expMinusNuX = std::exp( -eta * tanTheta );
      //       double temp = 1 + expMinusTwoX + 2. * b * expMinusX;

      //       double tempOne = std::sqrt(( expMinusX + temp * qK * qF * vNaught * vNaught ) / expMinusX);
      //       double common = std::sin(eta * PI) * (1. - expMinusTwoX) / temp
      //                        * (1 + tanTheta * tanTheta)
      //                        * kernel( tau, std::acosh(tempOne) ) / tempOne;
      //       if ( !(strike < EPSILON) )
      //         return common * std::cosh(eta * tanTheta);
      //       else
      //         return common * std::sinh(eta * tanTheta);
      //     };
      //     result += integrationMethod()->quadrature( fTwo, 0, .49 * PI );
      //     return factor * result + std::max(forward - strike, 0.);
      //   }
      //   opt_ptr_t updateModelParameters( const dbl_vec_t& params ) const override
      //   {
      //     return std::make_shared<FreeBoundarySABR>(params(0), params(1), params(2), params(3), kernelMethod(), integrationMethod());
      //   }
      // private:
      //   double doEffectiveStrike( double strike,
      //                             double forward )  const override
      //   {
      //     return std::max( strike, .1 * forward );
      //   }
      // };

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

    // beagle::pricer_ptr_t formClosedFormExactSABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
    //                                                                  const beagle::real_function_ptr_t& discounting,
    //                                                                  const beagle::real_function_ptr_t& alpha,
    //                                                                  const beagle::real_function_ptr_t& beta,
    //                                                                  const beagle::real_function_ptr_t& rho,
    //                                                                  const beagle::real_function_ptr_t& nu,
    //                                                                  bool useApproximateKernel,
    //                                                                  const beagle::integration_method_ptr_t& quadMethod)
    // {
    //   return std::make_shared<impl::ClosedFormExactSABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, useApproximateKernel, quadMethod);
    // }

    // beagle::pricer_ptr_t formClosedFormFreeBoundarySABREuropeanOptionPricer(const beagle::real_function_ptr_t& forward,
    //                                                                         const beagle::real_function_ptr_t& discounting,
    //                                                                         const beagle::real_function_ptr_t& alpha,
    //                                                                         const beagle::real_function_ptr_t& beta,
    //                                                                         const beagle::real_function_ptr_t& rho,
    //                                                                         const beagle::real_function_ptr_t& nu,
    //                                                                         bool useApproximateKernel,
    //                                                                         const beagle::integration_method_ptr_t& quadMethod)
    // {
    //   return std::make_shared<impl::ClosedFormFreeBoundarySABREuropeanOptionPricer>(forward, discounting, alpha, beta, rho, nu, useApproximateKernel, quadMethod);
    // }

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
