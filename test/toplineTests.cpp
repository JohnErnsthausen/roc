#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>
#include "data.hpp"
#include "matrix.hpp"
#include "vectorf.hpp"
#include "topline.hpp"

using namespace testing;

class TestThatTopLine : public Test
{
 public:
  const int num_coeff{30};
  double scale = 1.0;
  double rc = 0.0;
  double order = 0.0;
  std::vector<double> coeff;

  double epsilon{DBL_EPSILON};

  void SetUp() override
  {
    // Initialized to zero
    coeff = std::vector<double>(num_coeff);
  }

  void TearDown() override {}
};

TEST_F(TestThatTopLine, ThrowsExceptionIfToplineCalledWithCoeffSizeLessThanTOPLINE_KSTARTPlusTOPLINE_NUSE)
{
  coeff = std::vector<double>(TOPLINE_KSTART+TOPLINE_NUSE-1);
  ASSERT_THROW(topline(coeff, scale, rc, order), std::exception);
}

TEST_F(TestThatTopLine, ThrowsExceptionIfConstructLinearLeastSquaresSystemCalledWithWRowsUnequalToCoeffSizeMinusTOPLINE_KSTART)
{
  int ml{(int)coeff.size() - TOPLINE_KSTART - 1}, nl{2};
  matrix<double> WL(ml, nl);
  vectorf<double> betaL(ml);
  
  EXPECT_THROW(constructLinearLeastSquaresSystem(coeff, TOPLINE_KSTART, WL, betaL), std::exception);

  int mg{(int)coeff.size() - TOPLINE_KSTART + 1}, ng{2};
  matrix<double> WG(mg, ng);
  vectorf<double> betaG(mg);
  
  EXPECT_THROW(constructLinearLeastSquaresSystem(coeff, TOPLINE_KSTART, WG, betaG), std::exception);
}

TEST_F(TestThatTopLine, ThrowsExceptionIfLinearLeastSquaresSystemCalledWithUnequalWRowsAndBRows)
{
  int m{(int)coeff.size() - TOPLINE_KSTART}, n{4};
  matrix<double> W(m, n);
  vectorf<double> beta(m-1);

  ASSERT_THROW(constructLinearLeastSquaresSystem(coeff, TOPLINE_KSTART, W, beta), std::exception);
}

TEST_F(TestThatTopLine, ConstructsLinearLeastSquaresRow)
{
  double w1{0.0}, w2{0.0}, b{0.0};
  int k{3};
  coeff[3] = 1.0;

  constructLinearLeastSquaresRow(coeff, k, w1, w2, b);
  EXPECT_THAT(w1, DoubleEq(1.0));
  EXPECT_THAT(w2, DoubleEq(3.0));
  EXPECT_THAT(b, DoubleEq(0.0));
}

TEST_F(TestThatTopLine, CanSolveAContrivedLinearLeastSquaresSystem)
{
  for (int k{0}; k < num_coeff; k++)
  {
    coeff[ k ] = std::pow(10.0, -3.0 * (double)k + 1.0);
  }

  double w1, w2, b;
  for (int k{TOPLINE_KSTART}; k < (int)coeff.size(); k++)
  {
    constructLinearLeastSquaresRow(coeff, k, w1, w2, b);
    EXPECT_THAT(w1, DoubleNear(1.0, epsilon));
    EXPECT_THAT(w2, DoubleNear((double)k, epsilon));
    EXPECT_THAT(b, DoubleNear(-3.0 * (double)k + 1.0, epsilon));
  }
  
  double error = topline(coeff, scale, rc, order);
  EXPECT_THAT(error, DoubleNear(0.0, epsilon));

  // rc = 1/10^{-3} = 1000
  EXPECT_THAT(rc, DoubleNear(1000.0, 2.16006e-12));
}

// Problem from
// https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx
// Rc should be exactly 1/8.
//
// We want to fit the tail of the Taylor series, omit first
// TOPLINE_KSTART terms from consideration.
TEST_F(TestThatTopLine, WillComputeLeastSquaresSolutionOnExample)
{
  for (int k = 0; k < num_coeff; k++)
  {
    coeff[ k ] = std::pow(8, k + 1) / ((double)(k + 1));
  }

  EXPECT_THAT(topline(coeff, scale, rc, order),
              DoubleNear(0.0, 0.000966623));
  EXPECT_THAT(rc, DoubleNear(1.0/8.0, 0.131588));
  EXPECT_TRUE(std::isnan(order));
}
