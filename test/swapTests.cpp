#include <gmock/gmock.h>
#include <gtest/gtest.h>

extern "C"
{
#include "swap.h"
}

using namespace testing;

TEST(TestThatThisProject, HasSwapDoubleMethod)
{
  double a = 1.0, b = 2.0;

  dswap(&a, &b);

  EXPECT_THAT(a, DoubleEq(2.0));
  EXPECT_THAT(b, DoubleEq(1.0));
}

TEST(TestThatThisProject, HasSwapIntegerMethod)
{
  int a = 1, b = 2;

  iswap(&a, &b);

  EXPECT_THAT(a, Eq(2));
  EXPECT_THAT(b, Eq(1));
}
