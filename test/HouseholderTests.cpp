#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "utils.h"

#define DOUBLE_NEAR(x) DoubleNear((x), 4.0*DBL_EPSILON)
#define DOUBLE_NEAR_MULTIPLIER(x,multiplier) DoubleNear((x), (multiplier)*DBL_EPSILON)
#define map(i,j) ((j)-1)*lda+((i)-1)
#define a(i,j) a[map(i,j)]
#define r(i,j) r[map(i,j)]
#define map1(i) (i)-1
#define x(i) x[map1(i)]
#define y(i) y[map1(i)]
#define ipiv(i) ipiv[map1(i)]
#define wrk(i) wrk[map1(i)]
#define tau(i) tau[map1(i)]

using namespace std;
using namespace testing;

TEST(TestThatUtils, SgnFunctionReturnsNegativeOneWhenEvaluatedAtNegativeTwelve)
{
  ASSERT_THAT(sgn(-12.0), Eq(-1));
}

TEST(TestThatUtils, SgnFunctionReturnsPositiveOneWhenEvaluatedAtPositiveTwelve)
{
  ASSERT_THAT(sgn(12.0), Eq(1.0));
}

TEST(TestThatUtils, SgnFunctionReturnsZeroWhenEvaluatedAtZero)
{
  ASSERT_THAT(sgn(0.0), Eq(0.0));
}

TEST(TestThatUtils, SignFunctionReturnsTwelveWhenEvaluatedAtNegativeTwelveAndOne)
{
  ASSERT_THAT(sign(-12.0,1.0), DoubleEq(12.0));
}

TEST(TestThatUtils, SignFunctionReturnsTwelveWhenEvaluatedAtNegativeTwelveAndZero)
{
  ASSERT_THAT(sign(-12.0,0.0), DoubleEq(12.0));
}

TEST(TestThatUtils, SignFunctionReturnsNegativeTwelveWhenEvaluatedAtNegativeTwelveAndNegativeOne)
{
  ASSERT_THAT(sign(-12.0,-1.0), DoubleEq(-12.0));
}

class TestThatDDIST2: public Test
{
 public:
  int n=10, incx=1, incy=1, kvec=2;
  double *xvector = NULL;
  double *yvector = NULL;

  void SetUp() override
  {
    xvector = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (xvector == NULL)
    {
      printf("Failed to allocate xvector memory!\n");
      assert(xvector);
    }

    yvector = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (yvector == NULL)
    {
      printf("Failed to allocate yvector memory!\n");
      assert(yvector);
    }
  }

  void TearDown() override
  {
    free(yvector);
    free(xvector);
  }
};

TEST_F(TestThatDDIST2, ReturnsZeroWheneverTheDimensionNIsNotPositive)
{
  n= 0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
  n=-1;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
  n=-10;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2, ReturnsZeroWheneverKVECEqualsAnyIntegerExceptTwoAndINCXEqualsZero)
{
  kvec= 1;
  incx= 0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2, ReturnsZeroWheneverKVECEqualsTwoAndINCXEqualsZeroAndINCYEqualsZero)
{
  kvec= 2;
  incx= 0;
  incy= 0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2, ReturnsZeroWheneverXVectorIsZeroAndYVectorIsZeroAtDefaultInputs)
{
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2, ReturnsDistanceForSpecifiedXVectorAndYVector)
{
  xvector[2]= 3.0;
  yvector[2]= 2.999999999999981;

  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(1.898481372109018e-14));

  xvector[2]= 3.0;
  yvector[2]= 2.999999999999981;
  xvector[7]= 5.0;
  yvector[7]= 2.0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(3.0));
}

TEST_F(TestThatDDIST2, ReturnsDistanceForSpecifiedXVector)
{
  xvector[2]= 3.0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(3.0));

  xvector[2]= 3.0;
  xvector[7]= 3.0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(4.242640687119286));
}

TEST_F(TestThatDDIST2, HasNoUnderflowNoOverflow)
{
  const double cutlo = 4.44089e-16;
  const double cuthi = 1.30438e19;

  xvector[0]=1.0;
  xvector[1]=1.0;
  xvector[2]= cutlo;
  xvector[3]=1.0;
  xvector[4]=1.0;
  xvector[5]= cuthi;
  xvector[6]=1.0;
  xvector[7]= cutlo;
  xvector[8]=1.0;

  yvector[0]=1.0;
  yvector[1]=1.0;
  yvector[2]= cuthi;
  yvector[3]=1.0;
  yvector[4]=1.0;
  yvector[5]= cutlo;
  yvector[6]=1.0;
  yvector[7]= cuthi;
  yvector[8]=1.0;
  EXPECT_THAT(ddist2(n-1,xvector,incx,yvector,incy,kvec), DOUBLE_NEAR(2.259252432376692e+19));
}

class TestThatHouseholderG: public Test
{
 public:
  int n=10, incx=1;
  double alpha=2.0;
  double tau=2.0;
  double *xvector = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;
  const double cutlo = 4.44089e-16;
  const double cuthi = 1.30438e19;

  void SetUp() override
  {
    xvector = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (xvector == NULL)
    {
      printf("Failed to allocate xvector memory!\n");
      assert(xvector);
    }
  }

  void TearDown() override
  {
    free(xvector);
  }
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
  int myIncx=0;
  EXPECT_THAT(housg(n, &alpha, xvector, myIncx, &tau, safmin), Eq(-2));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholderG, IsDefinedForNegativeDimensionsButFlagsNegativeOne)
{
  int dim=-10;
  EXPECT_THAT(housg(dim, &alpha, xvector, incx, &tau, safmin), Eq(-1));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholderG, IsIdentityForOneDimensionalApplications)
{
  int dim=1;
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
  xvector[0]=1.0;
  xvector[1]=1.0;
  xvector[2]=1.0;
  xvector[3]=1.0;
  xvector[4]=1.0;
  xvector[5]=1.0;
  xvector[6]=1.0;
  xvector[7]=1.0;
  xvector[8]=1.0;
  EXPECT_THAT(ddist2(n-1,xvector,incx,xvector,incx,1), DOUBLE_NEAR(3.0));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatHouseholderG, ProductOfHTransposeHIsIdentityAndVerifyNormNoUnderflow)
{
  xvector[0]= cutlo/9;
  xvector[1]= cutlo/8;
  xvector[2]= cutlo/7;
  xvector[3]= cutlo/6;
  xvector[4]= cutlo/5;
  xvector[5]= cutlo/4;
  xvector[6]= cutlo/3;
  xvector[7]= cutlo/2;
  xvector[8]= cutlo;
  EXPECT_THAT(ddist2(n-1,xvector,incx,xvector,incx,1), DOUBLE_NEAR(5.510583948830441e-16));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatHouseholderG, ProductOfHTransposeHIsIdentityAndVerifyNormNoOverflow)
{
  xvector[0]= cuthi*2;
  xvector[1]= cuthi*3;
  xvector[2]= cuthi*4;
  xvector[3]= cuthi*5;
  xvector[4]= cuthi*6;
  xvector[5]= cuthi*7;
  xvector[6]= cuthi*8;
  xvector[7]= cuthi*9;
  xvector[8]= cuthi*10;
  EXPECT_THAT(ddist2(n-1,xvector,incx,xvector,incx,1), DOUBLE_NEAR(2.556052344553217e+20+32768));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatHouseholderG, ProductOfHTransposeHIsIdentityAndVerifyNormNoUnderflowNoOverflow)
{
  xvector[0]=1.0;
  xvector[1]=1.0;
  xvector[2]= cutlo;
  xvector[3]=1.0;
  xvector[4]=1.0;
  xvector[5]= cuthi;
  xvector[6]=1.0;
  xvector[7]= cutlo;
  xvector[8]=1.0;
  EXPECT_THAT(ddist2(n-1,xvector,incx,xvector,incx,1), DOUBLE_NEAR(1.304380000000000e+19));

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  double vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));
}

class TestThatHouseholderL: public Test
{
 public:
  int m=4, lda=4, n=3, ier=0, incx=1;
  double tau=2.0;
  double *a = NULL;
  double *r = NULL;
  double *xvector = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;

  void SetUp() override
  {
    xvector = (double *) calloc((size_t) m, (size_t) sizeof(double));
    if (xvector == NULL)
    {
      printf("Failed to allocate xvector memory!\n");
      assert(xvector);
    }

    a = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (a == NULL)
    {
      printf("Failed to allocate matrix A memory!\n");
      assert(a);
    }

    a(1,1) = 1.0, a(1,2) =  2.0, a(1,3) =  3.0; 
    a(2,1) = 1.0, a(2,2) =  5.0, a(2,3) =  6.0;
    a(3,1) = 1.0, a(3,2) =  8.0, a(3,3) =  9.0;
    a(4,1) = 1.0, a(4,2) = 11.0, a(4,3) = 12.0;

    r = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (r == NULL)
    {
      printf("Failed to allocate space for result r!\n");
      assert(r);
    }

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
  m=-1;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL, IsNotDefinedForNegativeColumnDimensionN)
{
  n=-1;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL, IsNotDefinedForINCXEqualsZero)
{
  incx=0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL, IsNotDefinedForLeadingDimensionLDALessThanNumberOfRowsM)
{
  lda=1;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatHouseholderL, IsDefinedForRowDimensionMEqualsZero)
{
  m=0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatHouseholderL, IsDefinedForColumnDimensionNEqualsZero)
{
  n=0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatHouseholderL, IsDefinedForTAUEqualsZero)
{
  tau=0.0;
  EXPECT_THAT(housl(m, n, xvector, incx, tau, a, lda, &ier), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatHouseholderL, ImplementsHouseholderMultiplicationOnLeft)
{
  xvector[0]= a(1,1);
  xvector[1]= a(2,1);
  xvector[2]= a(3,1);
  xvector[3]= a(4,1);

  double *alpha = xvector;
  double *v = xvector+1;

  // On input of housg, alpha is the first component of xvector
  EXPECT_THAT(*alpha, DoubleEq(xvector[0]));
  // On input of housg, v is components 2:m of xvector
  for(int row=0; row<m-1; row++)
    EXPECT_THAT(*(v+row), DoubleEq(xvector[row+1]));

  EXPECT_THAT(housg(m, alpha, v, incx, &tau, safmin), Eq(0));

  // On output of housg, beta is stored in alpha
  double beta = *alpha;
  // On output of housg, beta is stored in the first component of xvector
  EXPECT_THAT(beta, DoubleEq(xvector[0]));
  // On output of housg, xvector is components 2:m of v
  for(int row=0; row<m-1; row++)
    EXPECT_THAT(*(v+row), DoubleEq(xvector[row+1]));
  
  int i=1;
  EXPECT_THAT(housl(m-i+1, n-i, xvector, incx, tau, a+map(i,i+1), lda, &ier), Eq(0));

  // First column remains a(:,1) here which does not hold (beta v)
  r(1,1) = 1.0, r(1,2) =-13.0, r(1,3) =-15.0; 
  r(2,1) = 1.0, r(2,2) =  0.0, r(2,3) =  0.0;
  r(3,1) = 1.0, r(3,2) =  3.0, r(3,3) =  3.0;
  r(4,1) = 1.0, r(4,2) =  6.0, r(4,3) =  6.0;

  for(int col=1; col<=n; col++)
    for(int row=1; row<=m; row++)
      EXPECT_THAT(a(row,col), DOUBLE_NEAR_MULTIPLIER(r(row,col),100));
}

TEST_F(TestThatHouseholderL, ImplementsHouseholderMultiplicationOnLeftWithoutAdditionalStorage)
{
  double *alpha = a + map(1,1);
  double *v = a + map(2,1);

  // On input of housg, alpha is the first component of xvector
  EXPECT_THAT(*alpha, DoubleEq(a(1,1)));
  // On input of housg, v is components 2:m of xvector
  for(int row=0; row<m-1; row++)
    EXPECT_THAT(*(v+row), DoubleEq(a(row+2,1)));

  EXPECT_THAT(housg(m, alpha, v, incx, &tau, safmin), Eq(0));

  // On output of housg, beta is stored in alpha
  double beta = *alpha;
  // On output of housg, beta is stored in the first component of xvector
  EXPECT_THAT(beta, DoubleEq(a(1,1)));
  // On output of housg, xvector is components 2:m of v
  for(int row=0; row<m-1; row++)
    EXPECT_THAT(*(v+row), DoubleEq(a(row+2,1)));
  
  int i=1;
  EXPECT_THAT(housl(m-i+1, n-i, a+map(i,i), incx, tau, a+map(i,i+1), lda, &ier), Eq(0));

  // First column remains a(:,1) holds (beta v)
  r(1,1) =   beta, r(1,2) =-13.0, r(1,3) =-15.0; 
  r(2,1) = *(v+0), r(2,2) =  0.0, r(2,3) =  0.0;
  r(3,1) = *(v+1), r(3,2) =  3.0, r(3,3) =  3.0;
  r(4,1) = *(v+2), r(4,2) =  6.0, r(4,3) =  6.0;

  for(int col=1; col<=n; col++)
    for(int row=1; row<=m; row++)
      EXPECT_THAT(a(row,col), DOUBLE_NEAR_MULTIPLIER(r(row,col),100));
}

class TestThatQRF: public Test
{
 public:
  int m=4, lda=4, n=3, min_m_n=0, ier=0;
  double *a = NULL;
  int *ipiv = NULL;
  double *tau = NULL;
  double *wrk = NULL;
  double *r = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;

  void SetUp() override
  {
    a = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (a == NULL)
    {
      printf("Failed to allocate matrix A memory!\n");
      assert(a);
    }

    a(1,1) = 1.0, a(1,2) =  2.0, a(1,3) =  3.0; 
    a(2,1) = 1.0, a(2,2) =  5.0, a(2,3) =  6.0;
    a(3,1) = 1.0, a(3,2) =  8.0, a(3,3) =  9.0;
    a(4,1) = 1.0, a(4,2) = 11.0, a(4,3) = 12.0;

    min_m_n = (m<n) ? m : n;

    tau = (double *) calloc((size_t) min_m_n, (size_t) sizeof(double));
    if (tau == NULL)
    {
      printf("Failed to allocate tau memory!\n");
      assert(tau);
    }

    wrk = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (wrk == NULL)
    {
      printf("Failed to allocate wrk memory!\n");
      assert(wrk);
    }

    ipiv = (int *) calloc((size_t) n, (size_t) sizeof(int));
    if (ipiv == NULL)
    {
      printf("Failed to allocate ipiv memory!\n");
      assert(ipiv);
    }

    r = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (r == NULL)
    {
      printf("Failed to allocate space for result r!\n");
      assert(r);
    }
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
  ASSERT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(0));
}

TEST_F(TestThatQRF, IsNotDefinedForNegativeRowDimensionM)
{
  m=-1;
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRF, IsNotDefinedForNegativeColumnDimensionN)
{
  n=-1;
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRF, IsNotDefinedForLeadingDimensionLDALessThanNumberOfRowsM)
{
  lda=1;
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRF, IsDefinedForRowDimensionMEqualsZero)
{
  m=0;
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRF, IsDefinedForColumnDimensionNEqualsZero)
{
  n=0;
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRF, WillQRFactorA)
{
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(0));
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
//     -1.825741858350554e-01    -8.164965809277261e-01     5.477210117775342e-01
//     -3.651483716701108e-01    -4.082482904638630e-01    -7.312645785255080e-01
//     -5.477225575051662e-01    -8.814441016525298e-17    -1.806338782815866e-01
//     -7.302967433402215e-01     4.082482904638629e-01     3.641774450295606e-01
//   Column 4
//      1.301252240837732e-03
//      4.065121353587260e-01
//     -8.169280274399654e-01
//      4.091146398404016e-01
// R =
//     -1.643167672515498e+01    -1.825741858350554e+00    -1.460593486680443e+01
//                          0    -8.164965809277261e-01     8.164965809277228e-01
//                          0                         0     2.174446367394293e-15
//                          0                         0                         0
// E =
//      0     1     0
//      0     0     1
//      1     0     0
// >> transpose(Q)*Q
// ans =
//   Columns 1 through 3
//      1.000000000000000e+00     1.110223024625157e-16    -7.307522642552300e-17
//      1.110223024625157e-16     1.000000000000000e+00    -1.110223024625157e-16
//     -7.307522642552300e-17    -1.110223024625157e-16     1.000000000000000e+00
//                          0     8.283304597789254e-17     8.326672684688674e-17
//   Column 4
//                          0
//      8.283304597789254e-17
//      8.326672684688674e-17
//      9.999999999999999e-01
// >> Q*R-M*E
// ans =
//     -4.440892098500626e-16                         0     4.440892098500626e-15
//                          0                         0     8.881784197001252e-16
//                          0                         0     1.776356839400250e-15
//                          0                         0     1.776356839400250e-15
  // Check R 
  r(1,1) =-1.643167672515498e+01; r(1,2) =-1.825741858350554e+00; r(1,3) =-1.460593486680443e+01;
  r(2,1) =                     0; r(2,2) =-8.164965809277261e-01; r(2,3) = 8.164965809277228e-01;
  r(3,1) =                     0; r(3,2) =                     0; r(3,3) = 2.174446367394293e-15;
  r(4,1) =                     0; r(4,2) =                     0; r(4,3) =                     0;
  for(int col=1; col<=n; col++)
    for(int row=1; row<=col; row++)
      EXPECT_THAT(a(row,col), DOUBLE_NEAR_MULTIPLIER(r(row,col),100));
  // Check E
  EXPECT_THAT(ipiv[0], Eq(3));
  EXPECT_THAT(ipiv[1], Eq(1));
  EXPECT_THAT(ipiv[2], Eq(2));
}

class TestThatQRS: public Test
{
 public:
  int m=4, lda=4, n=2, min_m_n=0, ier=0;
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
    a = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (a == NULL)
    {
      printf("Failed to allocate matrix A memory!\n");
      assert(a);
    }

    a(1,1) = 1.0, a(1,2) =  2.0; 
    a(2,1) = 1.0, a(2,2) =  5.0;
    a(3,1) = 1.0, a(3,2) =  8.0;
    a(4,1) = 1.0, a(4,2) = 11.0;

    min_m_n = (m<n) ? m : n;

    tau = (double *) calloc((size_t) min_m_n, (size_t) sizeof(double));
    if (tau == NULL)
    {
      printf("Failed to allocate tau memory!\n");
      assert(tau);
    }

    wrk = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (wrk == NULL)
    {
      printf("Failed to allocate wrk memory!\n");
      assert(wrk);
    }

    ipiv = (int *) calloc((size_t) n, (size_t) sizeof(int));
    if (ipiv == NULL)
    {
      printf("Failed to allocate ipiv memory!\n");
      assert(ipiv);
    }

    r = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (r == NULL)
    {
      printf("Failed to allocate space for result r!\n");
      assert(r);
    }

    x = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (x == NULL)
    {
      printf("Failed to allocate space for result x!\n");
      assert(x);
    }

    y = (double *) calloc((size_t) m, (size_t) sizeof(double));
    if (y == NULL)
    {
      printf("Failed to allocate space for result y!\n");
      assert(y);
    }
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
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRS, IsNotDefinedForNegativeRowDimensionM)
{
  m=-1;
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsNotDefinedForNegativeColumnDimensionN)
{
  n=-1;
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsNotDefinedForRowDimensionMLessThanNumberColumnDimensionN)
{
  m=1;
  n=2;
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsNotDefinedForLeadingDimensionLDALessThanNumberOfRowsM)
{
  lda=1;
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(1));
  EXPECT_THAT(ier, Eq(1));
}

TEST_F(TestThatQRS, IsDefinedForRowDimensionMEqualsZero)
{
  m=0;
  n=0;
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRS, IsDefinedForColumnDimensionNEqualsZero)
{
  n=0;
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
}

TEST_F(TestThatQRS, WillQRFactorAANDComputeLinearLeastSquarsSolution)
{
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
  y(1)=3.0;
  y(2)=5.0;
  y(3)=2.0;
  y(4)=1.0;
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
  EXPECT_THAT(x(2), DOUBLE_NEAR_MULTIPLIER(4.7,100.0));
  EXPECT_THAT(x(1), DOUBLE_NEAR_MULTIPLIER(-0.3,100.0));
}

class TestThatQRSIntegrationTest: public Test
{
 public:
  int m=13, lda=13, n=2, min_m_n=0, ier=0;
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
    a = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (a == NULL)
    {
      printf("Failed to allocate matrix A memory!\n");
      assert(a);
    }

    a( 1,1) = 1.0, a( 1,2) = 18.0; 
    a( 2,1) = 1.0, a( 2,2) = 19.0;
    a( 3,1) = 1.0, a( 3,2) = 20.0;
    a( 4,1) = 1.0, a( 4,2) = 21.0;
    a( 5,1) = 1.0, a( 5,2) = 22.0; 
    a( 6,1) = 1.0, a( 6,2) = 23.0;
    a( 7,1) = 1.0, a( 7,2) = 24.0;
    a( 8,1) = 1.0, a( 8,2) = 25.0;
    a( 9,1) = 1.0, a( 9,2) = 26.0; 
    a(10,1) = 1.0, a(10,2) = 27.0;
    a(11,1) = 1.0, a(11,2) = 28.0;
    a(12,1) = 1.0, a(12,2) = 29.0;
    a(13,1) = 1.0, a(13,2) = 30.0; 

    min_m_n = (m<n) ? m : n;

    tau = (double *) calloc((size_t) min_m_n, (size_t) sizeof(double));
    if (tau == NULL)
    {
      printf("Failed to allocate tau memory!\n");
      assert(tau);
    }

    wrk = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (wrk == NULL)
    {
      printf("Failed to allocate wrk memory!\n");
      assert(wrk);
    }

    ipiv = (int *) calloc((size_t) n, (size_t) sizeof(int));
    if (ipiv == NULL)
    {
      printf("Failed to allocate ipiv memory!\n");
      assert(ipiv);
    }

    r = (double *) calloc((size_t) (lda*n), (size_t) sizeof(double));
    if (r == NULL)
    {
      printf("Failed to allocate space for result r!\n");
      assert(r);
    }

    x = (double *) calloc((size_t) n, (size_t) sizeof(double));
    if (x == NULL)
    {
      printf("Failed to allocate space for result x!\n");
      assert(x);
    }

    y = (double *) calloc((size_t) m, (size_t) sizeof(double));
    if (y == NULL)
    {
      printf("Failed to allocate space for result y!\n");
      assert(y);
    }

    y( 1) = log10(abs(pow(8,18)/18));
    y( 2) = log10(abs(pow(8,19)/19));
    y( 3) = log10(abs(pow(8,20)/20));
    y( 4) = log10(abs(pow(8,21)/21));
    y( 5) = log10(abs(pow(8,22)/22));
    y( 6) = log10(abs(pow(8,23)/23));
    y( 7) = log10(abs(pow(8,24)/24));
    y( 8) = log10(abs(pow(8,25)/25));
    y( 9) = log10(abs(pow(8,26)/26));
    y(10) = log10(abs(pow(8,27)/27));
    y(11) = log10(abs(pow(8,28)/28));
    y(12) = log10(abs(pow(8,29)/29));
    y(13) = log10(abs(pow(8,30)/30));
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

// Problem from
// https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx
//
// We want to fit the tail of the Taylor series, ignore the first 17 terms of the TCs.
//
// Solve on the first two components
// transpose(Q)*b=transpose(Q)*A*E*E*x=transpose(Q)*Q*R*E*x=R*E*x
TEST_F(TestThatQRSIntegrationTest, WillComputeLeastSquarsSolution)
{
  EXPECT_THAT(qrf( m, n, a, lda, ipiv, tau, wrk, safmin, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
  EXPECT_THAT(qrs( m, n, a, lda, tau, y, x, &ier ), Eq(0));
  EXPECT_THAT(ier, Eq(0));
  EXPECT_THAT(x(ipiv(1)), DOUBLE_NEAR_MULTIPLIER(-9.340358069665058e-01,100.0));
  EXPECT_THAT(x(ipiv(2)), DOUBLE_NEAR_MULTIPLIER( 8.847241983085994e-01,100.0));

  // If Taylor coefficients TC(i) are unscaled and T(i) = (1/scale)^i TC(i), then
  // multiply Rc by scale in this computation on scaled TCs.
  //
  // See this from ratio test as Rc is 1/L where
  //
  //   L = (1/scale)*limit of abs(TC(i+1)/TC(i)) = limit of abs(T(i+1)/T(i)).
  //
  // log10(abs(T(i))) = y(i) = b + m*i = [1 i](b,m) finds m and b.
  // Log10(L) = Log10(abs(T(i+1))/abs(T(i))) = Log10(abs(T(i+1))) - Log10(abs(T(i)))
  //           = m*(i+1)+b - (m*i+b) = m
  // rc = 1/L = scale/pow(10, m).
  double Rc = 1.0/pow(10,x(ipiv(2)));
  EXPECT_THAT(Rc, DOUBLE_NEAR_MULTIPLIER( 1.303994626296309e-01,100.0)); // Should be exactly 1/8
}

