#include "real_2d_function.hpp"
#include "real_function.hpp"

namespace beagle
{
  namespace math
  {
    namespace impl
    {
      struct TwoDimConstantFunction : public RealTwoDimFunction
      {
        TwoDimConstantFunction( double constant ) :
          m_Const( constant )
        { }
        virtual ~TwoDimConstantFunction( void )
        { }
      public:
        virtual double value( double argX,
                              double argY ) const override
        {
          return m_Const;
        }
      private:
        double m_Const;
      };

      struct PiecewiseConstantRightFunction : public RealTwoDimFunction
      {
        PiecewiseConstantRightFunction( const beagle::dbl_vec_t& params,
                                        const beagle::real_function_ptr_coll_t& funcs ) :
          m_Params( params ),
          m_Funcs( funcs )
        { }
        virtual ~PiecewiseConstantRightFunction( void )
        { }
      public:
        virtual double value( double argX,
                              double argY ) const override
        {
          if (argX >= m_Params.back())
            return m_Funcs.back()->value(argY);
          else
          {
            auto it = m_Params.cbegin();
            auto itEnd = m_Params.cend();
            auto jt = m_Funcs.cbegin();
            while (it != itEnd)
            {
              if (argX < *it)
                break;
              else
              {
                ++it;
                ++jt;
              }
            }

            return (*jt)->value(argY);
          }
        }
      private:
        beagle::dbl_vec_t m_Params;
        beagle::real_function_ptr_coll_t m_Funcs;
      };
    }


    RealTwoDimFunction::RealTwoDimFunction( void )
    { }

    RealTwoDimFunction::~RealTwoDimFunction( void )
    { }

    beagle::real_2d_function_ptr_t
    RealTwoDimFunction::createTwoDimConstantFunction( double constant )
    {
      return std::make_shared<impl::TwoDimConstantFunction>( constant );
    }

    beagle::real_2d_function_ptr_t
    RealTwoDimFunction::createPiecewiseConstantRightFunction( const beagle::dbl_vec_t& params,
                                                              const beagle::real_function_ptr_coll_t& funcs )
    {
      if (params.size() != funcs.size())
        throw(std::string("Cannot form piecewise constant function with unequal numbers of parameters and functions!"));

      return std::make_shared<impl::PiecewiseConstantRightFunction>( params, funcs );
    }
  }
}