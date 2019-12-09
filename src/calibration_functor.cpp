#include "calibration_functor.hpp"

#include <iostream>

namespace beagle
{
  namespace calibration
  {
    CalibrationFunctor::CalibrationFunctor( const dbl_vec_t& targets,
                                            const calibration_adapter_ptr_t& adapter,
                                            const Eigen::VectorXd& guesses,
                                            const calibration_bound_constraint_coll_t& constraints,
                                            const int_vec_t& calibIndices ) :
      m_Targets( targets ),
      m_Adapter( adapter ),
      m_Parameters( guesses ),
      m_Constraints( constraints ),
      m_CalibIndices( calibIndices )
    {
      if (guesses.size() != constraints.size())
        throw("Mismatch in the number of calibration guesses and calibration bound constraints!");
    }

    CalibrationFunctor::~CalibrationFunctor( void )
    { }

    int CalibrationFunctor::operator()( const Eigen::VectorXd& params, Eigen::VectorXd& diff ) const
    {
      dbl_vec_t parameters( params.size() );
      for (int i=0; i<parameters.size(); ++i)
      {
        parameters[i] = params(i);
      }

      dbl_vec_t values = m_Adapter->values(parameters, m_Constraints);
      for (int i=0; i<values.size(); ++i)
      {
        diff(i) = values[i] - m_Targets[i];
        std::cout << values[i] << "\t" << m_Targets[i] << "\t" << diff(i) << "\n";
      }

      std::cout << "\n\n";

      return 0;
    }
    int CalibrationFunctor::df(const Eigen::VectorXd& params, Eigen::MatrixXd& jacobian) const
    {
      dbl_vec_t parameters( params.size() );
      for (int i=0; i<parameters.size(); ++i)
      {
        parameters[i] = params(i);
      }

      dbl_mat_t derivs = m_Adapter->derivativeWithRespectToTransformedParameters( parameters, m_Constraints );

      for (int i=0; i<derivs.size(); ++i)
      {
        for (int j=0; j<m_Constraints.size(); ++j)
        {
          jacobian(i, j) = derivs[i][j];
        }
      }

      return 0;
    }
    int CalibrationFunctor::inputs( void ) const
    {
      return m_Parameters.size();
    }
    int CalibrationFunctor::values( void ) const
    {
      return m_Targets.size();
    }
    const Eigen::VectorXd& CalibrationFunctor::parameters( void ) const
    {
      return m_Parameters;
    }
  }
}
