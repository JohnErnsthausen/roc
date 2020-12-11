#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>
#include "topline.hpp"

using namespace testing;
using std::vector;

class TestThatTopLine : public Test
{
 public:
  const int num_tc{30};
  int kstart = 10;
  double scale = 1.0;
  double rc = 0.0, order = 0.0;
  vector<double> tc;
  double epsilon{DBL_EPSILON};

  void SetUp() override
  {
    tc = vector<double>(num_tc);
  }

  void TearDown() override { }
};

// Problem from
// https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx
// Rc should be exactly 1/8.
//
// We want to fit the tail of the Taylor series, ignore the first 17 terms of
// the TCs.
TEST_F(TestThatTopLine, WillComputeLeastSquarsSolution)
{
  for (int k = 0; k < num_tc; k++)
  {
    tc[ k ] = pow(8, k + 1) / ((double)(k + 1));
  }

  EXPECT_THAT(topline(tc, scale, rc, order), Eq(0));
  EXPECT_THAT(rc, DoubleNear(1.315873989428502e-01, epsilon));
}
