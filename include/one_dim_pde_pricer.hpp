#pragma once

#include "fwd_decl.hpp"
#include "pricer.hpp"

namespace beagle
{
  namespace valuation
  {
    struct OneDimParabolicPDEPricer : public Pricer
    {
      OneDimParabolicPDEPricer(const beagle::real_function_ptr_t& forward,
                               const beagle::real_function_ptr_t& discounting,
                               const beagle::real_2d_function_ptr_t& drift,
                               const beagle::real_2d_function_ptr_t& volatility,
                               const beagle::real_2d_function_ptr_t& rate,
                               const beagle::real_2d_function_ptr_t& recovery,
                               const beagle::valuation::OneDimFiniteDifferenceSettings& settings);
      virtual ~OneDimParabolicPDEPricer( void ) = default;
    protected:
      const beagle::real_function_ptr_t& forwardCurve( void ) const;
      const beagle::real_function_ptr_t& discountCurve( void ) const;
      const beagle::real_2d_function_ptr_t& volatilitySurface( void ) const;
      const beagle::valuation::OneDimFiniteDifferenceSettings& finiteDifferenceSettings(void) const;
      const beagle::real_2d_function_ptr_t& convectionCoefficient( void ) const;
      const beagle::real_2d_function_ptr_t& diffusionCoefficient( void ) const;
      const beagle::real_2d_function_ptr_t& rateCoefficient( void ) const;
      const beagle::real_2d_function_ptr_t& sourceTerm( void ) const;
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
