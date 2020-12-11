#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>
#include "data.hpp"
#include "matrix.hpp"
#include "topline.hpp"
#include "vectorf.hpp"

using namespace testing;
using std::vector;

class TestThatTopLine : public Test
{
 public:
  const int num_tc{30};
  double scale = 1.0;
  double rc = 0.0, order = 0.0;
  vector<double> tc;
  double epsilon{DBL_EPSILON};

  void SetUp() override { tc = vector<double>(num_tc); }

  void TearDown() override {}
};

TEST_F(TestThatTopLine, RequiresAtLeast10MoreThanTOPLINE_KSTARTCoefficients)
{
  tc = vector<double>(TOPLINE_KSTART + 9);
  ASSERT_THROW(topline(tc, scale, rc, order), std::exception);
}

TEST_F(TestThatTopLine, SatisfiesStorageRequirementsForLinearLeastSquaresSystem)
{
  int kstart{5};
  int m{num_tc - kstart - 1};
  int n{2};
  matrix<double> W(m, n);
  vectorf<double> b(m);

  ASSERT_THROW(constructLinearLeastSquaresSystem(tc, kstart, W, b),
               std::exception);
}

TEST_F(TestThatTopLine, CanConstructLinearLeastSquaresSystem)
{
  for (int k{0}; k < num_tc; k++)
  {
    tc[ k ] = (double)(k + 1);
  }

  int kstart{5};
  int m{num_tc - kstart};
  int n{2};
  matrix<double> W(m, n);
  vectorf<double> b(m);

  constructLinearLeastSquaresSystem(tc, kstart, W, b);

  for (int k{1}; k <= m; k++)
  {
    EXPECT_THAT(W(k, 1), DoubleNear(1.0, epsilon));
    EXPECT_THAT(W(k, 2), DoubleNear((double)(kstart + k - 1), epsilon));
    EXPECT_THAT(
        b(k), DoubleNear(log10(fabs((double)((kstart + k - 1) + 1))), epsilon));
  }
}

// Problem from
// https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx
// Rc should be exactly 1/8.
//
// We want to fit the tail of the Taylor series, omit first
// TOPLINE_KSTART terms from consideration.
TEST_F(TestThatTopLine, WillComputeLeastSquarsSolution)
{
  for (int k = 0; k < num_tc; k++)
  {
    tc[ k ] = pow(8, k + 1) / ((double)(k + 1));
  }

  EXPECT_THAT(topline(tc, scale, rc, order),
              DoubleNear(0.077610600696407323, epsilon));
  EXPECT_THAT(rc, DoubleNear(1.315873989428502e-01, epsilon));
  EXPECT_TRUE(std::isnan(order));
}
