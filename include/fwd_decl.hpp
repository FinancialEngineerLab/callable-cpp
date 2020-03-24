#ifndef FWD_DECL_HPP
#define FWD_DECL_HPP

#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <functional>

namespace beagle
{
  namespace product
  {
    struct Product;

    namespace option
    {
      struct Option;
      struct Payoff;
    }

    namespace bond
    {
      struct Bond;
    }
  }

  namespace math
  {
    struct RealFunction;
    struct RealTwoDimFunction;
    struct InterpolationBuilder;
    struct TwoDimInterpolationBuilder;
    struct OneDimParabolicPDESolver;
    struct BoundaryCondition;
  }

  namespace valuation
  {
    struct Pricer;
    struct DividendPolicy;
    struct FiniteDifferenceDetails;
  }

  namespace calibration
  {
    struct CalibrationBoundConstraint;
    struct CalibrationAdapter;
  }

  using product_ptr_t = std::shared_ptr<product::Product>;
  using product_ptr_coll_t = std::vector<product_ptr_t>;

  using payoff_ptr_t = std::shared_ptr<product::option::Payoff>;

  using real_function_ptr_t = std::shared_ptr<math::RealFunction>;
  using real_function_ptr_coll_t = std::vector<real_function_ptr_t>;
  using real_2d_function_ptr_t = std::shared_ptr<math::RealTwoDimFunction>;
  using interp_builder_ptr_t = std::shared_ptr<math::InterpolationBuilder>;
  using interp_builder_2d_ptr_t = std::shared_ptr<math::TwoDimInterpolationBuilder>;
  using parabolic_pde_solver_ptr_t = std::shared_ptr<math::OneDimParabolicPDESolver>;
  using boundary_condition_ptr_t = std::shared_ptr<math::BoundaryCondition>;

  using pricer_ptr_t = std::shared_ptr<valuation::Pricer>;
  using dividend_policy_ptr_t = std::shared_ptr<valuation::DividendPolicy>;

  using calibration_bound_constraint_ptr_t = std::shared_ptr<calibration::CalibrationBoundConstraint>;
  using calibration_bound_constraint_coll_t = std::vector<calibration_bound_constraint_ptr_t>;

  using calibration_adapter_ptr_t = std::shared_ptr<calibration::CalibrationAdapter>;

  using dbl_vec_t = std::vector<double>;
  using dbl_vec_vec_t = std::vector<dbl_vec_t>;
  using dbl_mat_t = dbl_vec_vec_t;
  using int_vec_t = std::vector<int>;

  using two_dbl_t = std::pair<double, double>;
  using discrete_dividend_schedule_t = std::vector<two_dbl_t>;

  using real_func_t = std::function<double(double)>;
  using real_2d_func_t = std::function<double(double, double)>;

  // payment time and cash amount with a face value of 100.0
  using bond_cashflows_t = std::vector<two_dbl_t>;
  using puttable_schedule_t = std::vector<two_dbl_t>;

  // For each callable schedule entry, the two doubles indicate the 
  // start and end of the call period.
  using callable_schedule_t = std::vector<std::tuple<double, beagle::real_function_ptr_t, double>>;
}



#endif
