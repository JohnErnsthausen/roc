#include <gmock/gmock.h>
#include <gtest/gtest.h>

extern "C"
{
#include "mathext.h"
}

using namespace testing;

TEST(TestThatMathext, SgnFunctionReturnsNegativeOneWhenEvaluatedAtNegativeTwelve)
{
  ASSERT_THAT(sgn(-12.0), Eq(-1));
}

TEST(TestThatMathext, SgnFunctionReturnsPositiveOneWhenEvaluatedAtPositiveTwelve)
{
  ASSERT_THAT(sgn(12.0), Eq(1.0));
}

TEST(TestThatMathext, SgnFunctionReturnsZeroWhenEvaluatedAtZero)
{
  ASSERT_THAT(sgn(0.0), Eq(0.0));
}

TEST(TestThatMathext,
     SignFunctionReturnsTwelveWhenEvaluatedAtNegativeTwelveAndOne)
{
  ASSERT_THAT(sign(-12.0, 1.0), DoubleEq(12.0));
}

TEST(TestThatMathext,
     SignFunctionReturnsTwelveWhenEvaluatedAtNegativeTwelveAndZero)
{
  ASSERT_THAT(sign(-12.0, 0.0), DoubleEq(12.0));
}

TEST(
    TestThatMathext,
    SignFunctionReturnsNegativeTwelveWhenEvaluatedAtNegativeTwelveAndNegativeOne)
{
  ASSERT_THAT(sign(-12.0, -1.0), DoubleEq(-12.0));
}
