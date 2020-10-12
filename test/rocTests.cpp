#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <cfloat>
#include <cassert>

extern "C"
{
#include "roc.h"
}

#define DOUBLE_NEAR(x) DoubleNear((x), 4.0 * DBL_EPSILON)
#define DOUBLE_NEAR_MULTIPLIER(x, multiplier) \
  DoubleNear((x), (multiplier)*DBL_EPSILON)

using namespace testing;

class TestThatROC : public Test
{
 public:
  int num_tc = 30, kstart = 18;
  double scale = 1.0;
  double rc = 0.0, slope=0.0, intercept=0.0;
  double *tc = NULL;

  void SetUp() override
  {
    tc = (double *)calloc((size_t)(num_tc), (size_t)sizeof(double));
    if (tc == NULL)
    {
      printf("Failed to allocate memory for Taylor coefficients!\n");
      assert(tc);
    }
  }

  void TearDown() override
  {
    free(tc);
  }
};

TEST_F(TestThatROC, CanAccessTheLeastSquaredSolutionMethod)
{
  EXPECT_THAT(roc(num_tc, tc, scale, kstart, &rc, &slope, &intercept), Eq(0));
}

// Problem from
// https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx
// Rc should be exactly 1/8.
//
// We want to fit the tail of the Taylor series, ignore the first 17 terms of
// the TCs.
TEST_F(TestThatROC, WillComputeLeastSquarsSolution)
{
  for(int k=0; k<num_tc; k++) { tc[k] = pow(8,k+1)/((double) (k+1)); }

  EXPECT_THAT(roc(num_tc, tc, scale, kstart, &rc, &slope, &intercept), Eq(0));
  EXPECT_THAT(intercept, DOUBLE_NEAR_MULTIPLIER(-9.340358069665058e-01, 100.0));
  EXPECT_THAT(slope, DOUBLE_NEAR_MULTIPLIER(8.847241983085994e-01, 100.0));
  EXPECT_THAT(rc, DOUBLE_NEAR_MULTIPLIER(1.303994626296309e-01, 100.0));
}
