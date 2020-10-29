#include "one_dim_pde_pricer.hpp"
#include "real_function.hpp"
#include "real_2d_function.hpp"

namespace beagle
{
  namespace valuation
  {
    OneDimParabolicPDEPricer::OneDimParabolicPDEPricer(const beagle::real_function_ptr_t& forward,
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
                        } );
    }

    const beagle::real_function_ptr_t& OneDimParabolicPDEPricer::forwardCurve( void ) const
    {
      return m_Forward;
    }

    const beagle::real_function_ptr_t& OneDimParabolicPDEPricer::discountCurve( void ) const
    {
      return m_Discounting;
    }

    const beagle::real_2d_function_ptr_t& OneDimParabolicPDEPricer::volatilitySurface( void ) const
    {
      return m_Vol;
    }

    const beagle::valuation::OneDimFiniteDifferenceSettings& OneDimParabolicPDEPricer::finiteDifferenceSettings(void) const
    {
      return m_Settings;
    }

    const beagle::real_2d_function_ptr_t& OneDimParabolicPDEPricer::convectionCoefficient( void ) const
    {
      return m_Convection;
    }

    const beagle::real_2d_function_ptr_t& OneDimParabolicPDEPricer::diffusionCoefficient( void ) const
    {
      return m_Diffusion;
    }

    const beagle::real_2d_function_ptr_t& OneDimParabolicPDEPricer::rateCoefficient( void ) const
    {
      return m_Rate;
    }

    const beagle::real_2d_function_ptr_t& OneDimParabolicPDEPricer::sourceTerm( void ) const
    {
      return m_Source;
    }
  }
}
