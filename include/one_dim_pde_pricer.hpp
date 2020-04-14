#ifndef ONE_DIM_PDE_PRICER_HPP
#define ONE_DIM_PDE_PRICER_HPP

#include "pricer.hpp"
#include "real_2d_function.hpp"
#include "option.hpp"
#include "payoff.hpp"
#include "util.hpp"
#include "dividend_policy.hpp"
#include "real_function.hpp"
#include "interpolation_builder.hpp"

namespace beagle
{
  namespace valuation
  {
    namespace impl
    {
      struct OneDimParabolicPDEPricer : public Pricer
      {
        OneDimParabolicPDEPricer(const beagle::real_function_ptr_t& forward,
                                 const beagle::real_function_ptr_t& discounting,
                                 const beagle::real_2d_function_ptr_t& drift,
                                 const beagle::real_2d_function_ptr_t& volatility,
                                 const beagle::real_2d_function_ptr_t& rate,
                                 const beagle::real_2d_function_ptr_t& recovery,
                                 const beagle::valuation::OneDimFiniteDifferenceSettings& settings) :
          m_Forward(forward),
          m_Discounting(discounting),
          m_Vol(volatility),
          m_Settings(settings)
        {
          m_Convection = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              double localVol = volatility->value(time, spot);
                              double dft = drift->value(time, spot);
                              return dft - .5 * localVol * localVol;
                            } );
          m_Diffusion = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              double localVol = volatility->value(time, spot);
                              return .5 * localVol * localVol;
                            } );
          m_Rate = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              return rate->value(time, spot);
                            } );
          m_Source = beagle::math::RealTwoDimFunction::createBinaryFunction(
                            [=](double time, double logMoneyness)
                            {
                              double spot = forward->value(time) * std::exp(logMoneyness);
                              return discounting->value(time) * recovery->value(time, spot);
                            } );}
        virtual ~OneDimParabolicPDEPricer( void ) = default;
      protected:
        const beagle::real_function_ptr_t& forwardCurve( void ) const
        {
          return m_Forward;
        }
        const beagle::real_function_ptr_t& discountCurve( void ) const
        {
          return m_Discounting;
        }
        const beagle::real_2d_function_ptr_t& volatilitySurface( void ) const
        {
          return m_Vol;
        }
        const beagle::valuation::OneDimFiniteDifferenceSettings& finiteDifferenceSettings(void) const
        {
          return m_Settings;
        }
        const beagle::real_2d_function_ptr_t& convectionCoefficient( void ) const
        {
          return m_Convection;
        }
        const beagle::real_2d_function_ptr_t& diffusionCoefficient( void ) const
        {
          return m_Diffusion;
        }
        const beagle::real_2d_function_ptr_t& rateCoefficient( void ) const
        {
          return m_Rate;
        }
        const beagle::real_2d_function_ptr_t& sourceTerm( void ) const
        {
          return m_Source;
        }
      private:
        beagle::real_function_ptr_t m_Forward;
        beagle::real_function_ptr_t m_Discounting;
        beagle::real_2d_function_ptr_t m_Vol;
        beagle::valuation::OneDimFiniteDifferenceSettings m_Settings;

        beagle::real_2d_function_ptr_t m_Convection;
        beagle::real_2d_function_ptr_t m_Diffusion;
        beagle::real_2d_function_ptr_t m_Rate;
        beagle::real_2d_function_ptr_t m_Source;
      };
    }
  }
}

#endif