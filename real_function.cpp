#include "real_function.hpp"

namespace beagle
{

  namespace impl
  {
    struct ConstantFunction : public RealFunction
    {
      ConstantFunction( double constant ) :
        m_Const( constant )
      { }
      virtual ~ConstantFunction( void )
      { }
    public:
      virtual double value( double arg ) const override
      {
        return m_Const;
      }
    private:
      double m_Const;
    };

    struct LinearWithFlatExtrapolationInterpolatedFunction : public RealFunction
    {
      LinearWithFlatExtrapolationInterpolatedFunction( const dbl_vec_t& xValues,
                                                       const dbl_vec_t& yValues ) :
        m_XValues( xValues ),
        m_YValues( yValues )
      { }
      virtual ~LinearWithFlatExtrapolationInterpolatedFunction( void )
      { }
    public:
      virtual double value( double arg ) const override
      {
        if (arg < m_XValues.front())
          return m_YValues.front();
        else if (arg >= m_XValues.back())
          return m_YValues.back();
        else
        {
          auto it = m_XValues.cbegin() + 1;
          auto itEnd = m_XValues.cend();
          auto jt = m_YValues.cbegin() + 1;
          while (it != itEnd)
          {
            if (arg < *it)
              break;
            else
            {
              ++it;
              ++jt;
            }
          }

          double xLeft = *(it - 1);
          double xRight = *it;
          double yLeft = *(jt - 1);
          double yRight = *jt;
          return yLeft + (arg - xLeft) / (xRight - xLeft) * (yRight - yLeft);
        }
      }
    private:
      dbl_vec_t m_XValues;
      dbl_vec_t m_YValues;
    };
  }


  RealFunction::RealFunction( void )
  { }

  RealFunction::~RealFunction( void )
  { }

  real_function_ptr_t
  RealFunction::createConstantFunction( double constant )
  {
    return std::make_shared<impl::ConstantFunction>( constant );
  }

  real_function_ptr_t
  RealFunction::createLinearWithFlatExtrapolationInterpolatedFunction( 
                                                       const dbl_vec_t& xValues,
                                                       const dbl_vec_t& yValues )
  {
    if (xValues.size() != yValues.size())
      throw(std::string("Mismatch in the size of interpolation parameters"));

    if (xValues.size() == 1U)
      return createConstantFunction( yValues[0] );
    else
      return std::make_shared<impl::LinearWithFlatExtrapolationInterpolatedFunction>( xValues,
                                                                                      yValues );
  }
}