#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <cassert>
#include <cfloat>

extern "C"
{
#include "dist.h"
}
#define DOUBLE_NEAR(x) DoubleNear((x), 4.0 * DBL_EPSILON)

using namespace testing;

class TestThatDDIST2 : public Test
{
 public:
  int n = 10, incx = 1, incy = 1, kvec = 2;
  double *xvector = NULL;
  double *yvector = NULL;

  void SetUp() override
  {
    xvector = (double *)calloc((size_t)n, (size_t)sizeof(double));
    if (xvector == NULL)
    {
      printf("Failed to allocate xvector memory!\n");
      assert(xvector);
    }

    yvector = (double *)calloc((size_t)n, (size_t)sizeof(double));
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
  n = 0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
  n = -1;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
  n = -10;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2,
       ReturnsZeroWheneverKVECEqualsAnyIntegerExceptTwoAndINCXEqualsZero)
{
  kvec = 1;
  incx = 0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2,
       ReturnsZeroWheneverKVECEqualsTwoAndINCXEqualsZeroAndINCYEqualsZero)
{
  kvec = 2;
  incx = 0;
  incy = 0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2,
       ReturnsZeroWheneverXVectorIsZeroAndYVectorIsZeroAtDefaultInputs)
{
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(0.0));
}

TEST_F(TestThatDDIST2, ReturnsDistanceForSpecifiedXVectorAndYVector)
{
  xvector[ 2 ] = 3.0;
  yvector[ 2 ] = 2.999999999999981;

  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec),
              DOUBLE_NEAR(1.898481372109018e-14));

  xvector[ 2 ] = 3.0;
  yvector[ 2 ] = 2.999999999999981;
  xvector[ 7 ] = 5.0;
  yvector[ 7 ] = 2.0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(3.0));
}

TEST_F(TestThatDDIST2, ReturnsDistanceForSpecifiedXVector)
{
  xvector[ 2 ] = 3.0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec), DOUBLE_NEAR(3.0));

  xvector[ 2 ] = 3.0;
  xvector[ 7 ] = 3.0;
  EXPECT_THAT(ddist2(n, xvector, incx, yvector, incy, kvec),
              DOUBLE_NEAR(4.242640687119286));
}

TEST_F(TestThatDDIST2, HasNoUnderflowNoOverflow)
{
  const double cutlo = 4.44089e-16;
  const double cuthi = 1.30438e19;

  xvector[ 0 ] = 1.0;
  xvector[ 1 ] = 1.0;
  xvector[ 2 ] = cutlo;
  xvector[ 3 ] = 1.0;
  xvector[ 4 ] = 1.0;
  xvector[ 5 ] = cuthi;
  xvector[ 6 ] = 1.0;
  xvector[ 7 ] = cutlo;
  xvector[ 8 ] = 1.0;

  yvector[ 0 ] = 1.0;
  yvector[ 1 ] = 1.0;
  yvector[ 2 ] = cuthi;
  yvector[ 3 ] = 1.0;
  yvector[ 4 ] = 1.0;
  yvector[ 5 ] = cutlo;
  yvector[ 6 ] = 1.0;
  yvector[ 7 ] = cuthi;
  yvector[ 8 ] = 1.0;
  EXPECT_THAT(ddist2(n - 1, xvector, incx, yvector, incy, kvec),
              DOUBLE_NEAR(2.259252432376692e+19));
}
