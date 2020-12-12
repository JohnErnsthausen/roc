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
