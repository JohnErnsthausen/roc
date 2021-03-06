#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <cassert>
#include <cfloat>

extern "C"
{
#include "dist.h"
#include "qrfactorization.h"
}

#define DOUBLE_NEAR(x) DoubleNear((x), 4.0 * DBL_EPSILON)
#define DOUBLE_NEAR_MULTIPLIER(x, multiplier) \
  DoubleNear((x), (multiplier)*DBL_EPSILON)
#define map(i, j) ((j)-1) * lda + ((i)-1)
#define a(i, j) a[ map(i, j) ]
#define q(i, j) q[ map(i, j) ]
#define r(i, j) r[ map(i, j) ]
#define map1(i) (i) - 1
#define x(i) x[ map1(i) ]
#define y(i) y[ map1(i) ]
#define ipiv(i) ipiv[ map1(i) ]
#define wrk(i) wrk[ map1(i) ]
#define tau(i) tau[ map1(i) ]

using namespace std;
using namespace testing;

class TestThatHouseholderG : public Test
{
 public:
  int n = 10, incx = 1;
  double alpha = 2.0;
  double tau = 2.0;
  double *xvector = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;
  const double cutlo = 4.44089e-16;
  const double cuthi = 1.30438e19;

  void SetUp() override
  {
    xvector = (double *)calloc((size_t)n, (size_t)sizeof(double));
    // if (xvector == NULL)
    //{
    //  printf("Failed to allocate xvector memory!\n");
    //  assert(xvector);
    //}
  }

  void TearDown() override { free(xvector); }
};

TEST_F(TestThatHouseholderG, CanAccessEPMACHandSAFMIN)
{
  EXPECT_THAT(0.2220446049250313E-15, DOUBLE_NEAR(epmach));
  EXPECT_THAT(0.2225073858507201E-307, DOUBLE_NEAR(safmin));
}

TEST_F(TestThatHouseholderG, CanAccessTheHouseholderReflectorMethod)
{
  ASSERT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
}

TEST_F(TestThatHouseholderG, IsDefinedForINCXEqualsZeroButFlagsNegativeTwo)
{
  int myIncx = 0;
  EXPECT_THAT(housg(n, &alpha, xvector, myIncx, &tau, safmin), Eq(-2));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholderG, IsDefinedForNegativeDimensionsButFlagsNegativeOne)
{
  int dim = -10;
  EXPECT_THAT(housg(dim, &alpha, xvector, incx, &tau, safmin), Eq(-1));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholderG, IsIdentityForOneDimensionalApplications)
{
  int dim = 1;
  EXPECT_THAT(housg(dim, &alpha, xvector, incx, &tau, safmin), Eq(0));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholderG, IsIdentityForXEqualsZero)
{
  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholderG, ProductOfHTransposeHIsIdentityAndVerifyNorm)
{
  xvector[ 0 ] = 1.0;
  xvector[ 1 ] = 1.0;
  xvector[ 2 ] = 1.0;
  xvector[ 3 ] = 1.0;
  xvector[ 4 ] = 1.0;
  xvector[ 5 ] = 1.0;
  xvector[ 6 ] = 1.0;
  xvector[ 7 ] = 1.0;
  xvector[ 8 ] = 1.0;
  EXPECT_THAT(ddist2(n - 1, xvector, incx, xvector, incx, 1), DOUBLE_NEAR(3.0));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2(n - 1, xvector, incx, xvector, incx, 1);
  EXPECT_THAT(-2.0 * tau + tau * tau * (1.0 + vnorm * vnorm), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatHouseholderG,
       ProductOfHTransposeHIsIdentityAndVerifyNormNoUnderflow)
{
  xvector[ 0 ] = cutlo / 9;
  xvector[ 1 ] = cutlo / 8;
  xvector[ 2 ] = cutlo / 7;
  xvector[ 3 ] = cutlo / 6;
  xvector[ 4 ] = cutlo / 5;
  xvector[ 5 ] = cutlo / 4;
  xvector[ 6 ] = cutlo / 3;
  xvector[ 7 ] = cutlo / 2;
  xvector[ 8 ] = cutlo;
  EXPECT_THAT(ddist2(n - 1, xvector, incx, xvector, incx, 1),
              DOUBLE_NEAR(5.510583948830441e-16));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2(n - 1, xvector, incx, xvector, incx, 1);
  EXPECT_THAT(-2.0 * tau + tau * tau * (1.0 + vnorm * vnorm), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatHouseholderG,
       ProductOfHTransposeHIsIdentityAndVerifyNormNoOverflow)
{
  xvector[ 0 ] = cuthi * 2;
  xvector[ 1 ] = cuthi * 3;
  xvector[ 2 ] = cuthi * 4;
  xvector[ 3 ] = cuthi * 5;
  xvector[ 4 ] = cuthi * 6;
  xvector[ 5 ] = cuthi * 7;
  xvector[ 6 ] = cuthi * 8;
  xvector[ 7 ] = cuthi * 9;
  xvector[ 8 ] = cuthi * 10;
  EXPECT_THAT(ddist2(n - 1, xvector, incx, xvector, incx, 1),
              DOUBLE_NEAR(2.556052344553217e+20 + 32768));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2(n - 1, xvector, incx, xvector, incx, 1);
  EXPECT_THAT(-2.0 * tau + tau * tau * (1.0 + vnorm * vnorm), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatHouseholderG,
       ProductOfHTransposeHIsIdentityAndVerifyNormNoUnderflowNoOverflow)
{
  xvector[ 0 ] = 1.0;
  xvector[ 1 ] = 1.0;
  xvector[ 2 ] = cutlo;
  xvector[ 3 ] = 1.0;
  xvector[ 4 ] = 1.0;
  xvector[ 5 ] = cuthi;
  xvector[ 6 ] = 1.0;
  xvector[ 7 ] = cutlo;
  xvector[ 8 ] = 1.0;
  EXPECT_THAT(ddist2(n - 1, xvector, incx, xvector, incx, 1),
              DOUBLE_NEAR(1.304380000000000e+19));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2(n - 1, xvector, incx, xvector, incx, 1);
  EXPECT_THAT(-2.0 * tau + tau * tau * (1.0 + vnorm * vnorm), DOUBLE_NEAR(0.0));
}

class TestThatHouseholderL : public Test
{
 public:
  int m = 4, lda = 4, n = 3, ier = 0, incx = 1;
  double tau = 2.0;
  double *a = NULL;
  double *r = NULL;
  double *xvector = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;

  void SetUp() override
  {
    xvector = (double *)calloc((size_t)m, (size_t)sizeof(double));
    // if (xvector == NULL)
    //{
    //  printf("Failed to allocate xvector memory!\n");
    //  assert(xvector);
    //}

    a = (double *)calloc((size_t)(lda * n), (size_t)sizeof(double));
    // if (a == NULL)
    //{
    //  printf("Failed to allocate matrix A memory!\n");
    //  assert(a);
    //}

    a(1, 1) = 1.0, a(1, 2) = 2.0, a(1, 3) = 3.0;
    a(2, 1) = 1.0, a(2, 2) = 5.0, a(2, 3) = 6.0;
    a(3, 1) = 1.0, a(3, 2) = 8.0, a(3, 3) = 9.0;
    a(4, 1) = 1.0, a(4, 2) = 11.0, a(4, 3) = 12.0;

    r = (double *)calloc((size_t)(lda * n), (size_t)sizeof(double));
    // if (r == NULL)
    //{
    //  printf("Failed to allocate space for result r!\n");
    //  assert(r);
    //}
  }

  void TearDown() override
  {
    free(r);
    free(a);
    free(xvector);
  }
};

TEST_F(TestThatHouseholderL, CanAccessTheHouseholderLeftMultiplyMethod)
{
  ASSERT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(0));
}

TEST_F(TestThatHouseholderL, IsNotDefinedForNegativeRowDimensionM)
{
  m = -1;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL, IsNotDefinedForNegativeColumnDimensionN)
{
  n = -1;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL, IsNotDefinedForINCXEqualsZero)
{
  incx = 0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL,
       IsNotDefinedForLeadingDimensionLDALessThanNumberOfRowsM)
{
  lda = 1;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL, IsDefinedForRowDimensionMEqualsZero)
{
  m = 0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatHouseholderL, IsDefinedForColumnDimensionNEqualsZero)
{
  n = 0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatHouseholderL, IsDefinedForTAUEqualsZero)
{
  tau = 0.0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatHouseholderL, ImplementsHouseholderMultiplicationOnLeft)
{
  xvector[ 0 ] = a(1, 1);
  xvector[ 1 ] = a(2, 1);
  xvector[ 2 ] = a(3, 1);
  xvector[ 3 ] = a(4, 1);

  double *alpha = xvector;
  double *v = xvector + 1;

  // On input of housg, alpha is the first component of xvector
  EXPECT_THAT(*alpha, DoubleEq(xvector[ 0 ]));
  // On input of housg, v is components 2:m of xvector
  for (int row = 0; row < m - 1; row++)
    EXPECT_THAT(*(v + row), DoubleEq(xvector[ row + 1 ]));

  EXPECT_THAT(housg(m, alpha, v, incx, &tau, safmin), Eq(0));

  // On output of housg, beta is stored in alpha
  double beta = *alpha;
  // On output of housg, beta is stored in the first component of xvector
  EXPECT_THAT(beta, DoubleEq(xvector[ 0 ]));
  // On output of housg, xvector is components 2:m of v
  for (int row = 0; row < m - 1; row++)
    EXPECT_THAT(*(v + row), DoubleEq(xvector[ row + 1 ]));

  int i = 1;
  EXPECT_THAT(
      housl(m - i + 1, n - i, xvector, incx, tau, a + map(i, i + 1), lda, &ier),
      Eq(0));

  // First column remains a(:,1) here which does not hold (beta v)
  r(1, 1) = 1.0, r(1, 2) = -13.0, r(1, 3) = -15.0;
  r(2, 1) = 1.0, r(2, 2) = 0.0, r(2, 3) = 0.0;
  r(3, 1) = 1.0, r(3, 2) = 3.0, r(3, 3) = 3.0;
  r(4, 1) = 1.0, r(4, 2) = 6.0, r(4, 3) = 6.0;

  for (int col = 1; col <= n; col++)
    for (int row = 1; row <= m; row++)
      EXPECT_THAT(a(row, col), DOUBLE_NEAR_MULTIPLIER(r(row, col), 100));
}

TEST_F(TestThatHouseholderL,
       ImplementsHouseholderMultiplicationOnLeftWithoutAdditionalStorage)
{
  double *alpha = a + map(1, 1);
  double *v = a + map(2, 1);

  // On input of housg, alpha is the first component of xvector
  EXPECT_THAT(*alpha, DoubleEq(a(1, 1)));
  // On input of housg, v is components 2:m of xvector
  for (int row = 0; row < m - 1; row++)
    EXPECT_THAT(*(v + row), DoubleEq(a(row + 2, 1)));

  EXPECT_THAT(housg(m, alpha, v, incx, &tau, safmin), Eq(0));

  // On output of housg, beta is stored in alpha
  double beta = *alpha;
  // On output of housg, beta is stored in the first component of xvector
  EXPECT_THAT(beta, DoubleEq(a(1, 1)));
  // On output of housg, xvector is components 2:m of v
  for (int row = 0; row < m - 1; row++)
    EXPECT_THAT(*(v + row), DoubleEq(a(row + 2, 1)));

  int i = 1;
  EXPECT_THAT(housl(m - i + 1, n - i, a + map(i, i), incx, tau,
                    a + map(i, i + 1), lda, &ier),
              Eq(0));

  // First column remains a(:,1) holds (beta v)
  r(1, 1) = beta, r(1, 2) = -13.0, r(1, 3) = -15.0;
  r(2, 1) = *(v + 0), r(2, 2) = 0.0, r(2, 3) = 0.0;
  r(3, 1) = *(v + 1), r(3, 2) = 3.0, r(3, 3) = 3.0;
  r(4, 1) = *(v + 2), r(4, 2) = 6.0, r(4, 3) = 6.0;

  for (int col = 1; col <= n; col++)
    for (int row = 1; row <= m; row++)
      EXPECT_THAT(a(row, col), DOUBLE_NEAR_MULTIPLIER(r(row, col), 100));
}

class TestThatQRF : public Test
{
 public:
  int m = 4, lda = 4, n = 3, min_m_n = 0, ier = 0;
  double *a = NULL;
  int *ipiv = NULL;
  double *tau = NULL;
  double *wrk = NULL;
  double *r = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;

  void SetUp() override
  {
    a = (double *)calloc((size_t)(lda * n), (size_t)sizeof(double));
    // if (a == NULL)
    //{
    //  printf("Failed to allocate matrix A memory!\n");
    //  assert(a);
    //}

    a(1, 1) = 1.0, a(1, 2) = 2.0, a(1, 3) = 3.0;
    a(2, 1) = 1.0, a(2, 2) = 5.0, a(2, 3) = 6.0;
    a(3, 1) = 1.0, a(3, 2) = 8.0, a(3, 3) = 9.0;
    a(4, 1) = 1.0, a(4, 2) = 11.0, a(4, 3) = 12.0;

    min_m_n = (m < n) ? m : n;

    tau = (double *)calloc((size_t)min_m_n, (size_t)sizeof(double));
    // if (tau == NULL)
    //{
    //  printf("Failed to allocate tau memory!\n");
    //  assert(tau);
    //}

    wrk = (double *)calloc((size_t)n, (size_t)sizeof(double));
    // if (wrk == NULL)
    //{
    //  printf("Failed to allocate wrk memory!\n");
    //  assert(wrk);
    //}

    ipiv = (int *)calloc((size_t)n, (size_t)sizeof(int));
    // if (ipiv == NULL)
    //{
    //  printf("Failed to allocate ipiv memory!\n");
    //  assert(ipiv);
    //}

    r = (double *)calloc((size_t)(lda * n), (size_t)sizeof(double));
    // if (r == NULL)
    //{
    //  printf("Failed to allocate space for result r!\n");
    //  assert(r);
    //}
  }

  void TearDown() override
  {
    free(r);
    free(ipiv);
    free(wrk);
    free(tau);
    free(a);
  }
};

TEST_F(TestThatQRF, CanAccessTheQRFactorizationMethod)
{
  ASSERT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(0));
}

TEST_F(TestThatQRF, IsNotDefinedForNegativeRowDimensionM)
{
  m = -1;
  EXPECT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRF, IsNotDefinedForNegativeColumnDimensionN)
{
  n = -1;
  EXPECT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRF, IsNotDefinedForLeadingDimensionLDALessThanNumberOfRowsM)
{
  lda = 1;
  EXPECT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRF, IsDefinedForRowDimensionMEqualsZero)
{
  m = 0;
  EXPECT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRF, IsDefinedForColumnDimensionNEqualsZero)
{
  n = 0;
  EXPECT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRF, WillQRFactorA)
{
  EXPECT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));

  // >> M = [[1 2 3];[1 5 6];[1 8 9];[1 11 12]]
  // M =
  //      1     2     3
  //      1     5     6
  //      1     8     9
  //      1    11    12
  // >> [Q,R,E] = qr(M)
  // Q =
  //   Columns 1 through 3
  //     -1.825741858350554e-01    -8.164965809277261e-01 5.477210117775342e-01
  //     -3.651483716701108e-01    -4.082482904638630e-01 -7.312645785255080e-01
  //     -5.477225575051662e-01    -8.814441016525298e-17 -1.806338782815866e-01
  //     -7.302967433402215e-01     4.082482904638629e-01 3.641774450295606e-01
  //   Column 4
  //      1.301252240837732e-03
  //      4.065121353587260e-01
  //     -8.169280274399654e-01
  //      4.091146398404016e-01
  // R =
  //     -1.643167672515498e+01    -1.825741858350554e+00 -1.460593486680443e+01
  //                          0    -8.164965809277261e-01 8.164965809277228e-01
  //                          0                         0 2.174446367394293e-15
  //                          0                         0 0
  // E =
  //      0     1     0
  //      0     0     1
  //      1     0     0
  // >> transpose(Q)*Q
  // ans =
  //   Columns 1 through 3
  //      1.000000000000000e+00     1.110223024625157e-16 -7.307522642552300e-17
  //      1.110223024625157e-16     1.000000000000000e+00 -1.110223024625157e-16
  //     -7.307522642552300e-17    -1.110223024625157e-16 1.000000000000000e+00
  //                          0     8.283304597789254e-17 8.326672684688674e-17
  //   Column 4
  //                          0
  //      8.283304597789254e-17
  //      8.326672684688674e-17
  //      9.999999999999999e-01
  // >> Q*R-M*E
  // ans =
  //     -4.440892098500626e-16                         0 4.440892098500626e-15
  //                          0                         0 8.881784197001252e-16
  //                          0                         0 1.776356839400250e-15
  //                          0                         0 1.776356839400250e-15
  // Check R
  r(1, 1) = -1.643167672515498e+01;
  r(1, 2) = -1.825741858350554e+00;
  r(1, 3) = -1.460593486680443e+01;
  r(2, 1) = 0;
  r(2, 2) = -8.164965809277261e-01;
  r(2, 3) = 8.164965809277228e-01;
  r(3, 1) = 0;
  r(3, 2) = 0;
  r(3, 3) = 2.174446367394293e-15;
  r(4, 1) = 0;
  r(4, 2) = 0;
  r(4, 3) = 0;
  for (int col = 1; col <= n; col++)
    for (int row = 1; row <= col; row++)
      EXPECT_THAT(a(row, col), DOUBLE_NEAR_MULTIPLIER(r(row, col), 100));
  // Check E
  EXPECT_THAT(ipiv[ 0 ], Eq(3));
  EXPECT_THAT(ipiv[ 1 ], Eq(1));
  EXPECT_THAT(ipiv[ 2 ], Eq(2));
}

class TestThatQRS : public Test
{
 public:
  int m = 4, lda = 4, n = 2, min_m_n = 0, ier = 0;
  double *a = NULL;
  int *ipiv = NULL;
  double *tau = NULL;
  double *wrk = NULL;
  double *r = NULL;
  double *x = NULL;
  double *y = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;

  void SetUp() override
  {
    a = (double *)calloc((size_t)(lda * n), (size_t)sizeof(double));
    // if (a == NULL)
    //{
    //  printf("Failed to allocate matrix A memory!\n");
    //  assert(a);
    //}

    a(1, 1) = 1.0, a(1, 2) = 2.0;
    a(2, 1) = 1.0, a(2, 2) = 5.0;
    a(3, 1) = 1.0, a(3, 2) = 8.0;
    a(4, 1) = 1.0, a(4, 2) = 11.0;

    min_m_n = (m < n) ? m : n;

    tau = (double *)calloc((size_t)min_m_n, (size_t)sizeof(double));
    // if (tau == NULL)
    //{
    //  printf("Failed to allocate tau memory!\n");
    //  assert(tau);
    //}

    wrk = (double *)calloc((size_t)n, (size_t)sizeof(double));
    // if (wrk == NULL)
    //{
    //  printf("Failed to allocate wrk memory!\n");
    //  assert(wrk);
    //}

    ipiv = (int *)calloc((size_t)n, (size_t)sizeof(int));
    // if (ipiv == NULL)
    //{
    //  printf("Failed to allocate ipiv memory!\n");
    //  assert(ipiv);
    //}

    r = (double *)calloc((size_t)(lda * n), (size_t)sizeof(double));
    // if (r == NULL)
    //{
    //  printf("Failed to allocate space for result r!\n");
    //  assert(r);
    //}

    x = (double *)calloc((size_t)n, (size_t)sizeof(double));
    // if (x == NULL)
    //{
    //  printf("Failed to allocate space for result x!\n");
    //  assert(x);
    //}

    y = (double *)calloc((size_t)m, (size_t)sizeof(double));
    // if (y == NULL)
    //{
    //  printf("Failed to allocate space for result y!\n");
    //  assert(y);
    //}
  }

  void TearDown() override
  {
    free(y);
    free(x);
    free(r);
    free(ipiv);
    free(wrk);
    free(tau);
    free(a);
  }
};

TEST_F(TestThatQRS, CanAccessTheLeastSquaredSolutionMethod)
{
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRS, IsNotDefinedForNegativeRowDimensionM)
{
  m = -1;
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsNotDefinedForNegativeColumnDimensionN)
{
  n = -1;
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsNotDefinedForRowDimensionMLessThanNumberColumnDimensionN)
{
  m = 1;
  n = 2;
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsNotDefinedForLeadingDimensionLDALessThanNumberOfRowsM)
{
  lda = 1;
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsDefinedForRowDimensionMEqualsZero)
{
  m = 0;
  n = 0;
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRS, IsDefinedForColumnDimensionNEqualsZero)
{
  n = 0;
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRS, WillQRFactorANDComputeLinearLeastSquarsSolution)
{
  EXPECT_THAT(qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
  y(1) = 3.0;
  y(2) = 5.0;
  y(3) = 2.0;
  y(4) = 1.0;
  EXPECT_THAT(qrs(m, n, a, lda, tau, y, x, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
  EXPECT_THAT(x(ipiv(1)), DOUBLE_NEAR_MULTIPLIER(4.7, 100.0));
  EXPECT_THAT(x(ipiv(2)), DOUBLE_NEAR_MULTIPLIER(-0.3, 100.0));
}

class TestThatQR : public Test
{
 public:
  int m{3}, lda{3}, n{3}, min_m_n{0}, ier{0}, ldq{3};
  vector<int> ipiv{0, 0, 0};
  vector<double> a{2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0};
  vector<double> q{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  vector<double> tau{0.0, 0.0, 0.0};
  vector<double> wrk{0.0, 0.0, 0.0};
  vector<double> x{0.0, 0.0, 0.0};
  vector<double> y{-1.0, -1.0, -1.0};
  double epmach{DBL_EPSILON};
  double safmin{DBL_MIN};

  vector<double> a_check{
      2.4494897427831783,  -0.57979589711327117, 0.28989794855663559,
      -1.6329931618554516, -1.5275252316519461,  -0.39985928207919952,
      -1.6329931618554518, 1.0910894511799623,   1.0690449676496974};
  vector<double> q_check{
      -0.40824829046386291, 0.81649658092772592,  -0.40824829046386296,
      -0.87287156094396945, -0.21821789023599258, 0.43643578047198495,
      0.26726124191242445,  0.53452248382484879,  0.80178372573727308};
  vector<int> ipiv_check{2, 1, 3};
  vector<double> x_check{-1.5, -2, -1.5};

  void SetUp() override { min_m_n = (m < n) ? m : n; }

  void TearDown() override {}
};

inline void EXPECT_DARRAY_EQ(const int n, double *expected, double *actual)
{
  for (int i = 0; i < n; i++)
  {
    EXPECT_THAT(*(actual + i), DOUBLE_NEAR_MULTIPLIER(*(expected + i), 100.0));
  }
}

inline void EXPECT_IARRAY_EQ(const int n, int *expected, int *actual)
{
  for (int i = 0; i < n; i++)
  {
    EXPECT_THAT(*(actual + i), Eq(*(expected + i)));
  }
}

TEST_F(TestThatQR, WillFactorANDComputeSolutionForSquareSystem)
{
  EXPECT_THAT(
      qrf(m, n, &a[ 0 ], lda, &ipiv[ 0 ], &tau[ 0 ], &wrk[ 0 ], safmin, &ier),
      Eq(0));
  EXPECT_THAT(ier, Eq(0));

  EXPECT_DARRAY_EQ(m * n, &a[ 0 ], &a_check[ 0 ]);
  EXPECT_IARRAY_EQ(m, &ipiv[ 0 ], &ipiv_check[ 0 ]);

  EXPECT_THAT(qorg(m, n, 1, m, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(0));
  EXPECT_THAT(ier, Eq(0));

  EXPECT_DARRAY_EQ(m * n, &a[ 0 ], &a_check[ 0 ]);
  EXPECT_DARRAY_EQ(m * n, &q[ 0 ], &q_check[ 0 ]);

  EXPECT_THAT(qrs(m, n, &a[ 0 ], lda, &tau[ 0 ], &y[ 0 ], &x[ 0 ], &ier),
              Eq(0));
  EXPECT_THAT(ier, Eq(0));

  EXPECT_DARRAY_EQ(m * n, &a[ 0 ], &a_check[ 0 ]);
  for (int i = 0; i < n; i++)
  {
    EXPECT_THAT(x[ ipiv[ i ] - 1 ],
                DOUBLE_NEAR_MULTIPLIER(x_check[ i ], 100.0));
  }
}

TEST_F(TestThatQR, WillProvideOrthogonalMartixQUnlessMLessThanZero)
{
  int ier{0};
  int m{-1};
  int n{2};
  int k{1};
  int ell{m};
  int lda{1};
  int ldq{2};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(1));
}

TEST_F(TestThatQR, WillProvideOrthogonalMartixQUnlessNLessThanZero)
{
  int ier{0};
  int m{2};
  int n{-1};
  int k{1};
  int ell{m};
  int lda{1};
  int ldq{2};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(1));
}

TEST_F(TestThatQR, WillProvideOrthogonalMartixQUnlessMLessThanN)
{
  int ier{0};
  int m{2};
  int n{3};
  int k{1};
  int ell{m};
  int lda{m};
  int ldq{m};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(1));
}

TEST_F(TestThatQR, WillProvideOrthogonalMartixQUnlessKGreaterThanL)
{
  int ier{0};
  int m{3};
  int n{2};
  int k{m + 1};
  int ell{m};
  int lda{m};
  int ldq{m};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(1));
}

TEST_F(TestThatQR, WillProvideOrthogonalMartixQUnlessLDALessThanM)
{
  int ier{0};
  int m{3};
  int n{2};
  int k{m};
  int ell{m};
  int lda{m - 1};
  int ldq{m};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(1));
}

TEST_F(TestThatQR, WillProvideOrthogonalMartixQUnlessLDQLessThanOne)
{
  int ier{0};
  int m{3};
  int n{2};
  int k{m};
  int ell{m};
  int lda{m};
  int ldq{0};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(1));
}

TEST_F(TestThatQR, WillProvideOrthogonalMartixQWithQuickReturnNEqualsZero)
{
  int ier{0};
  int m{3};
  int n{0};
  int k{m};
  int ell{m};
  int lda{m};
  int ldq{m};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(0));
}

TEST_F(TestThatQR,
       WillProvideOrthogonalMartixQWithQuickReturnNEqualsZeroANDMEqualsZero)
{
  int ier{0};
  int m{0};
  int n{0};
  int k{m};
  int ell{m};
  int lda{1};
  int ldq{1};
  // All matricies not accessed in this test
  EXPECT_THAT(qorg(m, n, k, ell, &a[ 0 ], lda, &tau[ 0 ], &q[ 0 ], ldq, &ier),
              Eq(0));
}
