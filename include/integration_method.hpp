#pragma once

#include "fwd_decl.hpp"

namespace beagle
{
  namespace math
  {
    struct IntegrationMethod
    {
      virtual ~IntegrationMethod(void) = default;
    public:
      virtual double quadrature(const real_function_ptr_t& f,
                                double lowerBound,
                                double upperBound) const = 0;
      virtual const int& numberOfSteps(void) const = 0;
    public:
      static integration_method_ptr_t midPointIntegrationMethod(int numSteps);
      static integration_method_ptr_t trapezoidIntegrationMethod(int numSteps);
    };
  }
}