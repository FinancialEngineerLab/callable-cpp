#include "gtest/gtest.h"
#include "util.hpp"
#include "interpolation_builder.hpp"
#include "real_function.hpp"

namespace beagle
{
  namespace test
  {
    TEST(test_tridiagonal_solver, Solve)
    {
      beagle::dbl_vec_t diag{ -2.6, -2.6, -2.6, -2.6 };
      beagle::dbl_vec_t upper{ 1., 1., 1., 1. };
      beagle::dbl_vec_t lower{ 1., 1., 1., 1. };
      beagle::dbl_vec_t rhs{ -240., 0., 0., -150. };

      beagle::dbl_vec_t results = {118.11216764581187, 67.09163587911088, 56.32608563987643, 79.35618678456785};
      beagle::util::tridiagonalSolve(rhs, diag, upper, lower);
      for (int i = 0; i < 4; ++i)
        EXPECT_NEAR(rhs[i], results[i], 1e-6);
    }

    TEST(test_natural_cubic_spline, Interpolation)
    {
      beagle::dbl_vec_t xValues{1., 2., 3., 4., 5.};
      beagle::dbl_vec_t yValues{0., 1., 0., 1., 0.};

      beagle::interp_builder_ptr_t spline = beagle::math::InterpolationBuilder::naturalCubicSpline();
      beagle::real_function_ptr_t func = spline->formFunction(xValues, yValues);

      for (int i = 0; i < xValues.size(); ++i)
        EXPECT_DOUBLE_EQ(func->value(xValues[i]), yValues[i]);

    }
  }
}