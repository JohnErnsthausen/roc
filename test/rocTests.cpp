#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>

#include "roc.hpp"

#define DOUBLE_NEAR(x) DoubleNear((x), DBL_EPSILON)

using namespace testing;
using namespace std;

class TestThatROC : public Test
{
};

TEST_F(TestThatROC, HasMethodDerivatives2Coeddicients)
{
  vector<double> coeffs{1,1,1,1,1};
  ASSERT_THAT(derivatives2coefficients(coeffs),Eq(0));
}

TEST_F(TestThatROC, MethodDerivatives2CoeddicientsConvertsCoefficientArrayOfDerivativesToCoefficientArrayOfCoefficients)
{
  vector<double> coeffs{1,1,1,1,1};
  derivatives2coefficients(coeffs);
  
  vector<double> expect{1,1,1.0/2,1.0/6,1.0/24};
  for(unsigned long int norder=1; norder<coeffs.size(); norder++)
  {
    EXPECT_THAT(coeffs[norder], DOUBLE_NEAR(expect[norder]));
  }
}

