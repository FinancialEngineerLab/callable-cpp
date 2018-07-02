#ifndef FWD_DECL_HPP
#define FWD_DECL_HPP

#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

namespace beagle
{
  namespace option
  {
    struct Option;
    struct Payoff;
  }

  namespace math
  {
    struct RealFunction;
    struct RealTwoDimFunction;
    struct InterpolationBuilder;
    struct TwoDimInterpolationBuilder;
  }

  namespace valuation
  {
    struct Pricer;
    struct DividendPolicy;
  }

  using option_ptr_t = std::shared_ptr<option::Option>;
  using payoff_ptr_t = std::shared_ptr<option::Payoff>;

  using real_function_ptr_t = std::shared_ptr<math::RealFunction>;
  using real_function_ptr_coll_t = std::vector<real_function_ptr_t>;
  using real_2d_function_ptr_t = std::shared_ptr<math::RealTwoDimFunction>;
  using interp_builder_ptr_t = std::shared_ptr<math::InterpolationBuilder>;
  using interp_builder_2d_ptr_t = std::shared_ptr<math::TwoDimInterpolationBuilder>;

  using pricer_ptr_t = std::shared_ptr<valuation::Pricer>;
  using dividend_policy_ptr_t = std::shared_ptr<valuation::DividendPolicy>;

  using dbl_vec_t = std::vector<double>;
  using int_vec_t = std::vector<int>;
  using discrete_dividend_schedule_t = std::vector< std::pair<double, double> >;
}



#endif