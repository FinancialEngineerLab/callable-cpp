#include "real_2d_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct BootstrappedLocalVolatilityFunction : public RealTwoDimFunction
      {
        BootstrappedLocalVolatilityFunction( const beagle::dbl_vec_t& expiries,
                                             const beagle::dbl_vec_vec_t& strikesColl,
                                             const beagle::dbl_vec_vec_t& pricesColl,
                                             const beagle::pricer_ptr_t& forwardPricer ) :
          m_Expiries( expiries ),
          m_StrikesColl( strikesColl ),
          m_PricesColl( pricesColl ),
          m_ForwardPricer( forwardPricer )
        { }
        virtual ~BootstrappedLocalVolatilityFunction( void )
        { }
      public:
        virtual double value( double argX,
                              double argY ) const override
        {
          if (!m_Func)
            doCalibration();

          return m_Func->value( argX, argY );
        }
      private:
        void doCalibration( void ) const
        {
          // TODO
        }
      private:
        beagle::dbl_vec_t m_Expiries;
        beagle::dbl_vec_vec_t m_StrikesColl;
        beagle::dbl_vec_vec_t m_PricesColl;
        beagle::pricer_ptr_t m_ForwardPricer;
        mutable beagle::real_2d_function_ptr_t m_Func;
      };

    }

    beagle::real_2d_function_ptr_t
    RealTwoDimFunction::createBootstrappedLocalVolatilityFunction( 
                                  const beagle::dbl_vec_t& expiries,
                                  const beagle::dbl_vec_vec_t& strikesColl,
                                  const beagle::dbl_vec_vec_t& pricesColl,
                                  const beagle::pricer_ptr_t& forwardPricer )
    {
      return std::make_shared<impl::BootstrappedLocalVolatilityFunction>( expiries,
                                                                          strikesColl,
                                                                          pricesColl,
                                                                          forwardPricer );
    }
  }
}