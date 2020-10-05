#pragma once

#include "fwd_decl.hpp"
#include "calibration_adapter.hpp"
#include "calibration_constraint.hpp"
#include "Eigen/Dense"

namespace beagle
{
  namespace calibration
  {
    template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
    struct Functor
    {
      typedef _Scalar Scalar;
      enum {
          InputsAtCompileTime = NX,
          ValuesAtCompileTime = NY
      };
      typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
      typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
      typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

      int m_inputs, m_values;

      Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
      Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

    };

    struct CalibrationFunctor : public Functor<double>
    {
      CalibrationFunctor( const dbl_vec_t& targets,
                          const calibration_adapter_ptr_t& adapter,
                          const Eigen::VectorXd& guesses,
                          const calibration_bound_constraint_coll_t& constraints,
                          const int_vec_t& calibIndices );
      virtual ~CalibrationFunctor( void ) = default;
    public:
      int operator()( const Eigen::VectorXd& params, Eigen::VectorXd& diff ) const;
      int df(const Eigen::VectorXd& params, Eigen::MatrixXd& jacobian) const;
      int inputs( void ) const;
      int values( void ) const;
      const Eigen::VectorXd& parameters( void ) const;
    private:
      dbl_vec_t m_Targets;
      calibration_adapter_ptr_t m_Adapter;
      Eigen::VectorXd m_Parameters;
      calibration_bound_constraint_coll_t m_Constraints;
      int_vec_t m_CalibIndices;
    };
  }
}
