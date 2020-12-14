#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <cfloat>
#include <vector>

#include "qrfactorization.hpp"

using namespace std;
using namespace testing;

class TestThatCPPQR : public Test
{
 public:
  int m{3}, n{3};
  double epsilon{DBL_EPSILON};

  void SetUp() override
  {

  }

  void TearDown() override {}
};

TEST_F(TestThatCPPQR, CanAccessTheQRFactorizationMethod)
{
  matrix<double> a(m, n);
  vectorf<double> x(n);
  vectorf<double> b(m);
  a(1, 1) = 2.0, a(1, 2) = -1.0, a(1, 3) = 0.0;
  a(2, 1) = -1.0, a(2, 2) = 2.0, a(2, 3) = -1.0;
  a(3, 1) = 0.0, a(3, 2) = -1.0, a(3, 3) = 2.0;
  b(1) = -1.0;
  b(2) = -1.0;
  b(3) = -1.0;
  ASSERT_THAT(qr(m, n, a, b, x), Eq(0));
}

TEST_F(TestThatCPPQR, WillThrowQRFactorizationErrorFromQRF)
{
  matrix<double> a(m, n, vector<double>(m*n, 0.0));
  vectorf<double> x(n);
  vectorf<double> b(m);
  // This throw comes from qrf finding a zero pivot
  //
  // Once there is a zero pivot found in qrf, the method qr
  // throws an error and returns.
  //
  // qrs will fail for the same reason. A zero pivot is the
  // only reason for a fail in qrs. The QRSolver error
  // throw is redundant. Code coverate should discover this.
  ASSERT_THROW(qr(m, n, a, b, x), std::exception);
}

TEST_F(TestThatCPPQR, WillReturnACorrectSolution)
{
  matrix<double> a(m, n);
  vectorf<double> x(n);
  vectorf<double> b(m);
  a(1, 1) = 2.0, a(1, 2) = -1.0, a(1, 3) = 0.0;
  a(2, 1) = -1.0, a(2, 2) = 2.0, a(2, 3) = -1.0;
  a(3, 1) = 0.0, a(3, 2) = -1.0, a(3, 3) = 2.0;
  b(1) = -1.0;
  b(2) = -1.0;
  b(3) = -1.0;
  qr(m, n, a, b, x);
  EXPECT_THAT(x(1), DoubleNear(-1.5, 8.88179e-16));
  EXPECT_THAT(x(2), DoubleNear(-2.0, epsilon));
  EXPECT_THAT(x(3), DoubleNear(-1.5, 8.88178e-16));
}

TEST_F(TestThatCPPQR, QRCanResolveAOneDimensionalSystem)
{
  matrix<double> a(1, 1);
  vectorf<double> x(1);
  vectorf<double> b(1);
  a(1, 1) = 2.0;
  b(1) = -1.0;

  qr(1, 1, a, b, x);
  
  ASSERT_THAT(x(1), DoubleNear(-0.5, epsilon));
}

TEST_F(TestThatCPPQR, QRFailsOnSingularOneDimensionalSystem)
{
  matrix<double> a(1, 1);
  vectorf<double> x(1);
  vectorf<double> b(1);
  a(1, 1) = 0.0;
  b(1) = -1.0;

  ASSERT_THROW(qr(1, 1, a, b, x), std::exception);
}

TEST_F(TestThatCPPQR, FactorFailsOnSingularOneDimensionalSystem)
{
  matrix<double> a(1, 1, vector<double> {0.0});
  vectorf<double> tau(1);
  vectorf<int> ipiv(1);

  ASSERT_THROW(factor(1, 1, a, tau, ipiv), std::exception);
}

TEST_F(TestThatCPPQR, FactorFailsOnSingularMultiDimensionalSystem)
{
  matrix<double> a(m, n);
  vectorf<double> tau(n);
  vectorf<int> ipiv(n);
  a(1, 1) = 2.0, a(1, 2) = -1.0, a(1, 3) = 0.0;
  a(2, 1) = -1.0, a(2, 2) = 2.0, a(2, 3) = -1.0;
  a(3, 1) = 0.0, a(3, 2) = 0.0, a(3, 3) = 0.0;

  ASSERT_THROW(factor(m, n, a, tau, ipiv), std::exception);
}

TEST_F(TestThatCPPQR, SolveFailsOnSingularOneDimensionalSystem)
{
  matrix<double> a(1, 1, vector<double> {0.0});
  vectorf<double> tau(1);
  vectorf<double> b(1);
  vectorf<double> x(1);

  ASSERT_THROW(solve(1, 1, a, tau, b, x), std::exception);
}

TEST_F(TestThatCPPQR, SolveFailsOnSingularMultiDimensionalSystem)
{
  matrix<double> a(m, n);
  vectorf<double> tau(n);
  vectorf<double> b(m);
  vectorf<double> x(n);
  a(1, 1) = 1.0, a(1, 2) = 0.0, a(1, 3) = 0.0;
  a(2, 1) = 0.0, a(2, 2) = 1.0, a(2, 3) = 0.0;
  a(3, 1) = 0.0, a(3, 2) = 0.0, a(3, 3) = 0.0;

  ASSERT_THROW(solve(m, n, a, tau, b, x), std::exception);
}

TEST_F(TestThatCPPQR, PermutesASolutionNVector)
{
  vectorf<int> ipiv(3);
  vectorf<double> wrk(3);
  vectorf<double> x(3);
  ipiv(1) = 3;
  ipiv(2) = 1;
  ipiv(3) = 2;
  x(1) = 2.0;
  x(2) = 3.0;
  x(3) = 1.0;
  
  permute(3,x,ipiv,wrk);
  EXPECT_THAT(x(1), DoubleNear(1.0, epsilon));
  EXPECT_THAT(x(2), DoubleNear(2.0, epsilon));
  EXPECT_THAT(x(3), DoubleNear(3.0, epsilon));
}




