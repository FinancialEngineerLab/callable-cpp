#include "integration_method.hpp"
#include "real_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct MidPointIntegrationMethod : public IntegrationMethod
      {
        MidPointIntegrationMethod(int numSteps) :
          m_NumSteps( numSteps )
        { }
      public:
        double quadrature( const real_function_ptr_t& f,
                           double lowerBound,
                           double upperBound ) const override
        {
          double step = (upperBound - lowerBound) / m_NumSteps;
          double result = 0.;

          for (int i = 0; i < m_NumSteps; ++i)
          {
            result += f->value(lowerBound + (i + .5) * step) * step;
          }

          return result;
        }
        const int& numberOfSteps( void ) const override
        {
          return m_NumSteps;
        }
      private:
        int m_NumSteps;
      };

      struct TrapezoidIntegrationMethod : public IntegrationMethod
      {
        TrapezoidIntegrationMethod( int numSteps) :
          m_NumSteps( numSteps )
        { }
      public:
        double quadrature( const real_function_ptr_t& f,
                           double lowerBound,
                           double upperBound ) const override
        {
          double step = (upperBound - lowerBound) / m_NumSteps;
          double arg = lowerBound;
          double result = .5 * f->value(arg);

          for (int i = 1; i < m_NumSteps; ++i)
          {
            arg += step;
            result += f->value(arg) * step;
          }

          return result + .5 * f->value(upperBound);
        }
        const int& numberOfSteps( void ) const override
        {
          return m_NumSteps;
        }
      private:
        int m_NumSteps;
      };
    }

    integration_method_ptr_t IntegrationMethod::midPointIntegrationMethod( int numSteps )
    {
      return std::make_shared<impl::MidPointIntegrationMethod>( numSteps );
    }

    integration_method_ptr_t IntegrationMethod::trapezoidIntegrationMethod( int numSteps )
    {
      return std::make_shared<impl::TrapezoidIntegrationMethod>( numSteps );
    }
  }

}