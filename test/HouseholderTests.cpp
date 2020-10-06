#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "utils.h"

#define DOUBLE_NEAR(x) DoubleNear((x), 1e-15)

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

TEST(TestThatUtils, HasSwapIntegerMethod)
{
  double a=1.0, b=2.0;

  swap(&a,&b);

  EXPECT_THAT(a,DoubleEq(2.0));
  EXPECT_THAT(b,DoubleEq(1.0));
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

class TestThatHouseholder: public Test
{
 public:
  int n=10, incx=1;
  double alpha=2.0;
  double tau=2.0;
  double *xvector = NULL;
  double epmach = DBL_EPSILON;
  double safmin = DBL_MIN;

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

TEST_F(TestThatHouseholder, CanAccessEPMACHandSAFMIN)
{
  EXPECT_THAT(0.2220446049250313E-15, DOUBLE_NEAR(epmach));
  EXPECT_THAT(0.2225073858507201E-307, DOUBLE_NEAR(safmin));
}

TEST_F(TestThatHouseholder, CanAccessTheHouseholderReflectorMethod)
{
  ASSERT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
}

TEST_F(TestThatHouseholder, IsDefinedForINCXEqualsZeroButFlagsNegativeTwo)
{
  int myIncx=0;
  EXPECT_THAT(housg(n, &alpha, xvector, myIncx, &tau, safmin), Eq(-2));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholder, IsDefinedForNegativeDimensionsButFlagsNegativeOne)
{
  int dim=-10;
  EXPECT_THAT(housg(dim, &alpha, xvector, incx, &tau, safmin), Eq(-1));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholder, IsIdentityForOneDimensionalApplications)
{
  int dim=1;
  EXPECT_THAT(housg(dim, &alpha, xvector, incx, &tau, safmin), Eq(0));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholder, IsIdentityForXEqualsZero)
{
  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  EXPECT_THAT(tau, DoubleEq(0.0));
}

TEST_F(TestThatHouseholder, ProductOfHTransposeHIsIdentity)
{
  const double cutlo = 4.44089e-16/n;
  const double cuthi = 1.30438e19*10.0;
  double vnorm = 0.0;

  // xvector components not less than cutlo or greater than cuthi in DDIST2
  xvector[0]=1.0;
  xvector[1]=1.0;
  xvector[2]=1.0;
  xvector[3]=1.0;
  xvector[4]=1.0;
  xvector[5]=1.0;
  xvector[6]=1.0;
  xvector[7]=1.0;
  xvector[8]=1.0;

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));

  xvector[0]= cutlo;
  xvector[1]= cutlo;
  xvector[2]= cutlo;
  xvector[3]= cutlo;
  xvector[4]= cutlo;
  xvector[5]= cutlo;
  xvector[6]= cutlo;
  xvector[7]= cutlo;
  xvector[8]= cutlo;

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));

  xvector[0]= cuthi;
  xvector[1]= cuthi;
  xvector[2]= cuthi;
  xvector[3]= cuthi;
  xvector[4]= cuthi;
  xvector[5]= cuthi;
  xvector[6]= cuthi;
  xvector[7]= cuthi;
  xvector[8]= cuthi;

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));

  xvector[0]=1.0;
  xvector[1]=1.0;
  xvector[2]= cutlo;
  xvector[3]=1.0;
  xvector[4]=1.0;
  xvector[5]= cuthi;
  xvector[6]=1.0;
  xvector[7]= cutlo;
  xvector[8]=1.0;

  EXPECT_THAT(housg(n, &alpha, xvector, incx, &tau, safmin), Eq(0));
  vnorm = ddist2( n-1,xvector,incx,xvector,incx,1);
  EXPECT_THAT(-2.0*tau + tau*tau*(1.0+vnorm*vnorm), DOUBLE_NEAR(0.0));
}

