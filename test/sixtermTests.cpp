#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>
#include <cmath>
#include <vector>
#include <iostream>

extern "C"
{
#include "qrfactorization.h"
}
#include "io.hpp"

#define TOL 1.0e-10
#define DOUBLE_NEAR(x) DoubleNear((x), TOL)

using namespace testing;
using namespace std;

// The six-term-test of Chang and Corliss
//
// The vector coeff is expected to have at least length 10
//
// The method returns ier=0 if computation was successful. Otherwise
// the method returns ier=1 whenever the algorithm detects that
// the coefficients do not resemble a pole.
int sixTerm(const vector<double> &coeff, const double &scale, double &rc,
            double &order)
{
  int m{4}, n{4}, ier{0};
  vector<int> ipiv(m, 0);
  vector<double> W(m*n, 0.0);
  vector<double> tau(m, 0.0);
  vector<double> wrk(m, 0.0);
  vector<double> x(m, 0.0);
  vector<double> b(m, 0.0);
  double safmin{DBL_MIN};

#define map1(i) (i)-1
#define b(i) b[ map1(i) ]
#define x(i) x[ map1(i) ]
#define ipiv(i) ipiv[ map1(i) ]
#define coeff(i) coeff[ map1(i) ]
#define map(i, j) ((j)-1) * m + ((i)-1)
#define W(i, j) W[ map(i, j) ]
  int nUse = 4;
  int k = coeff.size() - nUse;
  for(int i{1}; i<=nUse; i++)
  {
    b(i) = k * coeff(k+1);
    W(i,1) = 2.0 * coeff(k);
    W(i,2) = 2.0 * (k - 1) * coeff(k);
    W(i,3) =-2.0 * coeff(k-1);
    W(i,4) =-(k - 2) * coeff(k-1);
    k++;
  }

  qrf(m, n, &W[0], m, &ipiv[0], &tau[0], &wrk[0], safmin, &ier);
  if ( ier != 0 ) printf( "Solver error go to top-line\n" );
  qrs(m, n, &W[0], m, &tau[0], &b[0], &x[0], &ier);
  if ( ier != 0 ) printf( "Solver error go to top-line\n" );
  for(int i = 1; i <=n; i++) { b(i) = x(ipiv(i)); }

  double hOverRc{0.0}, cosTheta{0.0}, singularityOrder1{0.0}, singularityOrder2{0.0};

  // cout.precision(16);
  // cout << scientific;
  // cout << " bet12 = " << b(2) << '\n';
  // cout << " bet14 = " << b(4) << '\n';

  if( b(4) < 0 )
  {
    printf( "Runtime error go to top-line: Sqrt of negative number\n" );
    ier = 1;
    return ier;
  }

  // TODO Should the case b(4) == 0 be separated out to mean Rc is Inf?
  
  hOverRc = sqrt(b(4));
  rc  = scale / hOverRc;
  cosTheta = b(2)/hOverRc;

  // cout << " rc       = " << rc << '\n';
  // cout << " cosTheta = " << cosTheta << '\n';

  // Check -1 <= cosTheta <= 1
  if( (cosTheta < -1.0) || (cosTheta > 1.0) )
  {
    printf( "Runtime error go to top-line: CosTheta range\n" );
    ier = 1;
    return ier;
  }

  // // s is zero case
  // if( b(2) == 0 && b(4) == 0)
  // {
  //   order = 0.0;
  // }

  // // CosTheta can be zero
  // if( b(2) == 0 && b(4) != 0)
  // {
  //   singularityOrder2 = b(3)/b(4);
  //   order = singularityOrder2;
  // }

  // // cosTheta not zero and h/Rc not zero and b(4) is zero (Kind of strange)
  // if( b(2) != 0 && b(4) == 0)
  // {
  //   if( b(3) == 0 )
  //   {
  //     singularityOrder2 = b(3)/b(4);
  //     order = singularityOrder2;
  //   }
  //   else
  //   {
  //     printf( "Runtime error go to top-line\n" );
  //     ier = 1;
  //     return ier;
  //   }
  // }

  singularityOrder1 = b(1)/b(2);
  singularityOrder2 = b(3)/b(4);

  // cout << " s1 = " << singularityOrder1 << '\n';
  // cout << " s2 = " << singularityOrder2 << '\n';

  if( isnan(singularityOrder1) && isnan(singularityOrder2) )
  {
     printf( "Runtime error go to top-line: Order\n" );
     ier = 1;
     return ier;
  }
  
  if( isnan(singularityOrder1) && !isnan(singularityOrder2) )
    order = singularityOrder2;

  if( !isnan(singularityOrder1) && isnan(singularityOrder2) )
    order = singularityOrder1;

  if( !isnan(singularityOrder1) && !isnan(singularityOrder2) )
    order = (singularityOrder1+singularityOrder2)/2.0;

  // cout << " order = " << order << '\n';

  // Compare order
  // printf( "Abs of difference between two computations for order: [%22.16f]\n",
  //     fabs(singularityOrder1-singularityOrder2));
  
  // Check for agreement between computations against TOL at
  // previous W equation
  //
  // If no agreement (in backward error), check term analysis is said to have failed
  // becuase the coefficients do not represent a complex conjugate pair of poles.
  nUse++;
  k = coeff.size() - nUse;
  double check = k * coeff(k+1) - (
        (2.0 * coeff(k))*b(1) +
        (2.0 * (k - 1) * coeff(k))*b(2) +
        (-2.0 * coeff(k-1))*b(3) +
        (-(k - 2) * coeff(k-1))*b(4));
  if (fabs(check) > TOL)
  {
    cout.precision(16);
    cout << scientific;
    cout << "|check equation| = " << fabs(check) << "\n";
    ier = 1;
    return ier;
  }

  return ier;
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithScalingTenthForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t =-0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0/5.0};
  double time{-0.3};
  double scale{0.1};
  vector<double> coeffs{
    3.0769230769230771e-01,
    1.4201183431952663e-01,
    4.1875284478834776e-02,
    8.4030671194986195e-03,
    6.5716294139668780e-04,
    -3.4308380547065328e-04,
    -2.0889736724773902e-04,
    -7.0023107539675438e-05,
    -1.6249329076177969e-05,
    -2.1132974551840289e-06,
    2.7458033423644587e-07,
    2.8929072773866954e-07,
    1.1239723324581319e-07,
    2.9622513210477662e-08,
    5.0259881551579064e-09,
    4.1031978497675339e-11,
    -3.6767663724398882e-10,
    -1.7285321553550831e-10,
    -5.1495588920697002e-11,
    -1.0470793691436441e-11,
    -8.7147486368628090e-13,
    4.0322650071682737e-13,
    2.5314106676824962e-13,
    8.5816915376359253e-14,
    2.0135417345377379e-14,
    2.6919683612234631e-15,
    -3.0643132138743080e-16,
    -3.4850432996523446e-16,
    -1.3727651218492121e-16,
    -3.6550364857253289e-17 
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(0));
  EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  EXPECT_THAT(order, DOUBLE_NEAR(1.0));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithScalingTenthForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0/5.0};
  double time{0.0};
  double scale{0.1};
  vector<double> coeffs{
    1.0000000000000000e+00,
    0.0000000000000000e+00,
    -2.5000000000000000e-01,
    0.0000000000000000e+00,
    6.2500000000000000e-02,
    0.0000000000000000e+00,
    -1.5625000000000000e-02,
    0.0000000000000000e+00,
    3.9062500000000000e-03,
    0.0000000000000000e+00,
    -9.7656250000000000e-04,
    0.0000000000000000e+00,
    2.4414062500000000e-04,
    0.0000000000000000e+00,
    -6.1035156250000000e-05,
    0.0000000000000000e+00,
    1.5258789062500000e-05,
    0.0000000000000000e+00,
    -3.8146972656250000e-06,
    0.0000000000000000e+00,
    9.5367431640625000e-07,
    0.0000000000000000e+00,
    -2.3841857910156250e-07,
    0.0000000000000000e+00,
    5.9604644775390625e-08,
    0.0000000000000000e+00,
    -1.4901161193847656e-08,
    0.0000000000000000e+00,
    3.7252902984619141e-09,
    0.0000000000000000e+00 
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(0));
  EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  EXPECT_THAT(order, DOUBLE_NEAR(1.0));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithScalingTenthForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0/5.0};
  double time{0.3};
  double scale{0.1};
  vector<double> coeffs{
    3.0769230769230771e-01,
    -1.4201183431952663e-01,
    4.1875284478834776e-02,
    -8.4030671194986195e-03,
    6.5716294139668780e-04,
    3.4308380547065328e-04,
    -2.0889736724773902e-04,
    7.0023107539675438e-05,
    -1.6249329076177969e-05,
    2.1132974551840289e-06,
    2.7458033423644587e-07,
    -2.8929072773866954e-07,
    1.1239723324581319e-07,
    -2.9622513210477662e-08,
    5.0259881551579064e-09,
    -4.1031978497675339e-11,
    -3.6767663724398882e-10,
    1.7285321553550831e-10,
    -5.1495588920697002e-11,
    1.0470793691436441e-11,
    -8.7147486368628090e-13,
    -4.0322650071682737e-13,
    2.5314106676824962e-13,
    -8.5816915376359253e-14,
    2.0135417345377379e-14,
    -2.6919683612234631e-15,
    -3.0643132138743080e-16,
    3.4850432996523446e-16,
    -1.3727651218492121e-16,
    3.6550364857253289e-17 
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(0));
  EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  EXPECT_THAT(order, DOUBLE_NEAR(1.0));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t =-0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  //double a{1.0/5.0};
  //double time{-0.3};
  double scale{1.0};
  vector<double> coeffs{
   3.0769230769230771e-01,
    1.4201183431952664e+00,
    4.1875284478834782e+00,
    8.4030671194986173e+00,
    6.5716294139668623e+00,
    -3.4308380547065390e+01,
    -2.0889736724773917e+02,
    -7.0023107539675470e+02,
    -1.6249329076177976e+03,
    -2.1132974551840298e+03,
    2.7458033423644597e+03,
    2.8929072773866963e+04,
    1.1239723324581322e+05,
    2.9622513210477668e+05,
    5.0259881551579072e+05,
    4.1031978497674834e+04,
    -3.6767663724398911e+06,
    -1.7285321553550843e+07,
    -5.1495588920697041e+07,
    -1.0470793691436444e+08,
    -8.7147486368627921e+07,
    4.0322650071682835e+08,
    2.5314106676824994e+09,
    8.5816915376359320e+09,
    2.0135417345377384e+10,
    2.6919683612234612e+10,
    -3.0643132138743221e+10,
    -3.4850432996523499e+11,
    -1.3727651218492134e+12,
    -3.6550364857253306e+12
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(1));
  //EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  //EXPECT_THAT(order, DOUBLE_NEAR(1.0));
  //EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0+25.0*t*t);

  //double a{1.0/5.0};
  //double time{0.0};
  double scale{1.0};
  vector<double> coeffs{
    1.0000000000000000e+00,
    0.0000000000000000e+00,
    -2.5000000000000000e+01,
    0.0000000000000000e+00,
    6.2500000000000000e+02,
    0.0000000000000000e+00,
    -1.5625000000000000e+04,
    0.0000000000000000e+00,
    3.9062500000000000e+05,
    0.0000000000000000e+00,
    -9.7656250000000000e+06,
    0.0000000000000000e+00,
    2.4414062500000000e+08,
    0.0000000000000000e+00,
    -6.1035156250000000e+09,
    0.0000000000000000e+00,
    1.5258789062500000e+11,
    0.0000000000000000e+00,
    -3.8146972656250000e+12,
    0.0000000000000000e+00,
    9.5367431640625000e+13,
    0.0000000000000000e+00,
    -2.3841857910156250e+15,
    0.0000000000000000e+00,
    5.9604644775390624e+16,
    0.0000000000000000e+00,
    -1.4901161193847657e+18,
    0.0000000000000000e+00,
    3.7252902984619139e+19,
    0.0000000000000000e+00
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(1));
  //EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  //EXPECT_THAT(order, DOUBLE_NEAR(1.0));
  //EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  //double a{1.0/5.0};
  //double time{0.3};
  double scale{1.0};
  vector<double> coeffs{
    3.0769230769230771e-01,
    -1.4201183431952664e+00,
    4.1875284478834782e+00,
    -8.4030671194986173e+00,
    6.5716294139668623e+00,
    3.4308380547065390e+01,
    -2.0889736724773917e+02,
    7.0023107539675470e+02,
    -1.6249329076177976e+03,
    2.1132974551840298e+03,
    2.7458033423644597e+03,
    -2.8929072773866963e+04,
    1.1239723324581322e+05,
    -2.9622513210477668e+05,
    5.0259881551579072e+05,
    -4.1031978497674834e+04,
    -3.6767663724398911e+06,
    1.7285321553550843e+07,
    -5.1495588920697041e+07,
    1.0470793691436444e+08,
    -8.7147486368627921e+07,
    -4.0322650071682835e+08,
    2.5314106676824994e+09,
    -8.5816915376359320e+09,
    2.0135417345377384e+10,
    -2.6919683612234612e+10,
    -3.0643132138743221e+10,
    3.4850432996523499e+11,
    -1.3727651218492134e+12,
    3.6550364857253306e+12
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(1));
  //EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  //EXPECT_THAT(order, DOUBLE_NEAR(1.0));
  //EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithScalingHundredthForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t =-0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.14159};

  double a{1.0/5.0};
  double time{-0.3};
  double scale{0.01};
  vector<double> coeffs{
    2.4655761562767709e-02,
    3.5748957668969889e-03,
    2.8208219200667809e-04,
    1.5636812703033154e-05,
    6.5875429164410868e-07,
    2.1094322169896426e-08,
    4.5261475074311576e-10,
    1.1260909441922327e-12,
    -4.6867806776404976e-13,
    -2.8056792096691283e-14,
    -1.0573048948187840e-15,
    -2.8313453349783412e-17,
    -4.3638374553989396e-19,
    5.4965174709350763e-21,
    7.3086407830397857e-22,
    3.3112669205272192e-23,
    1.0201294809137287e-24,
    2.1125410354692833e-26,
    1.1958555283498346e-28,
    -1.3772058760031110e-29,
    -8.1538224313874953e-31,
    -2.8716158576770647e-32,
    -7.0504814826739304e-34,
    -9.3677299909588449e-36,
    1.6819506294660599e-37,
    1.6868283550256514e-38,
    6.9196677516458404e-40,
    1.9436081834894383e-41,
    3.5195594052044087e-43,
    2.8474916689612392e-46
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(0));
  EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  EXPECT_THAT(order, DOUBLE_NEAR(3.1415E0));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithScalingThousandthForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.14159};

  double a{1.0/5.0};
  double time{0.0};
  double scale{0.001};
  vector<double> coeffs{
    1.0000000000000000e+00,
    0.0000000000000000e+00,
    -7.8537500000000005e-05,
    0.0000000000000000e+00,
    4.0657882031250008e-09,
    0.0000000000000000e+00,
    -1.7420208371972664e-13,
    0.0000000000000000e+00,
    6.6866381072793826e-18,
    0.0000000000000000e+00,
    -2.3876313021567852e-22,
    0.0000000000000000e+00,
    8.0995417693789439e-27,
    0.0000000000000000e+00,
    -2.6443557530277716e-31,
    0.0000000000000000e+00,
    8.3805418341659806e-36,
    0.0000000000000000e+00,
    -2.5936613012600057e-40,
    0.0000000000000000e+00,
    7.8727346723120797e-45,
    0.0000000000000000e+00,
    -2.3513532430952043e-49,
    0.0000000000000000e+00,
    6.9274295598397387e-54,
    0.0000000000000000e+00,
    -2.0171475900060203e-58,
    0.0000000000000000e+00,
    5.8142478257289360e-63,
    0.0000000000000000e+00
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(0));
  EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  EXPECT_THAT(order, DOUBLE_NEAR(3.1415E0));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithScalingThousandthForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  double a{1.0/5.0};
  double time{0.3};
  double scale{0.001};
  vector<double> coeffs{
    2.4655761562767709e-02,
    -3.5748957668969888e-04,
    2.8208219200667818e-06,
    -1.5636812703033177e-08,
    6.5875429164411143e-11,
    -2.1094322169896626e-13,
    4.5261475074312502e-16,
    -1.1260909441924017e-19,
    -4.6867806776406748e-21,
    2.8056792096693676e-23,
    -1.0573048948189600e-25,
    2.8313453349792976e-28,
    -4.3638374554030017e-31,
    -5.4965174709220135e-34,
    7.3086407830372582e-36,
    -3.3112669205275378e-38,
    1.0201294809143237e-40,
    -2.1125410354727341e-43,
    1.1958555283637769e-46,
    1.3772058759989342e-48,
    -8.1538224313796743e-51,
    2.8716158576775504e-53,
    -7.0504814826864239e-56,
    9.3677299910270096e-59,
    1.6819506294408271e-61,
    -1.6868283550188946e-63,
    6.9196677516356937e-66,
    -1.9436081834912157e-68,
    3.5195594052256001e-71,
    -2.8474916699531912e-75
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(0));
  EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  EXPECT_THAT(order, DOUBLE_NEAR(3.1415E0));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  //double a{1.0/5.0};
  //double time{-0.3};
  double scale{1.0};
  vector<double> coeffs{
    2.4655761562767709e-02,
    3.5748957668969888e-01,
    2.8208219200667819e+00,
    1.5636812703033181e+01,
    6.5875429164411216e+01,
    2.1094322169896708e+02,
    4.5261475074313199e+02,
    1.1260909441928604e+02,
    -4.6867806776404268e+03,
    -2.8056792096692519e+04,
    -1.0573048948189148e+05,
    -2.8313453349791735e+05,
    -4.3638374554030015e+05,
    5.4965174709191988e+05,
    7.3086407830349784e+06,
    3.3112669205263112e+07,
    1.0201294809138143e+08,
    2.1125410354710954e+08,
    1.1958555283602265e+08,
    -1.3772058759989672e+09,
    -8.1538224313751421e+09,
    -2.8716158576746910e+10,
    -7.0504814826746704e+10,
    -9.3677299909914948e+10,
    1.6819506294477527e+11,
    1.6868283550187634e+12,
    6.9196677516269053e+12,
    1.9436081834863168e+13,
    3.5195594052075742e+13,
    2.8474916694776353e+12
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(1));
  //EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  //EXPECT_THAT(order, DOUBLE_NEAR(3.1415E0));
  //EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  //double a{1.0/5.0};
  //double time{0.0};
  double scale{1.0};
  vector<double> coeffs{
    1.0000000000000000e+00,
    0.0000000000000000e+00,
    -7.8537500000000009e+01,
    0.0000000000000000e+00,
    4.0657882031250006e+03,
    0.0000000000000000e+00,
    -1.7420208371972659e+05,
    0.0000000000000000e+00,
    6.6866381072793789e+06,
    0.0000000000000000e+00,
    -2.3876313021567836e+08,
    0.0000000000000000e+00,
    8.0995417693789396e+09,
    0.0000000000000000e+00,
    -2.6443557530277716e+11,
    0.0000000000000000e+00,
    8.3805418341659883e+12,
    0.0000000000000000e+00,
    -2.5936613012600116e+14,
    0.0000000000000000e+00,
    7.8727346723121130e+15,
    0.0000000000000000e+00,
    -2.3513532430952214e+17,
    0.0000000000000000e+00,
    6.9274295598398177e+18,
    0.0000000000000000e+00,
    -2.0171475900060536e+20,
    0.0000000000000000e+00,
    5.8142478257290687e+21,
    0.0000000000000000e+00
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(1));
  //EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  //EXPECT_THAT(order, DOUBLE_NEAR(3.1415E0));
  //EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  //double a{1.0/5.0};
  //double time{0.3};
  double scale{1.0};
  vector<double> coeffs{
    2.4655761562767709e-02,
    -3.5748957668969888e-01,
    2.8208219200667819e+00,
    -1.5636812703033181e+01,
    6.5875429164411216e+01,
    -2.1094322169896708e+02,
    4.5261475074313199e+02,
    -1.1260909441928604e+02,
    -4.6867806776404268e+03,
    2.8056792096692519e+04,
    -1.0573048948189148e+05,
    2.8313453349791735e+05,
    -4.3638374554030015e+05,
    -5.4965174709191988e+05,
    7.3086407830349784e+06,
    -3.3112669205263112e+07,
    1.0201294809138143e+08,
    -2.1125410354710954e+08,
    1.1958555283602265e+08,
    1.3772058759989672e+09,
    -8.1538224313751421e+09,
    2.8716158576746910e+10,
    -7.0504814826746704e+10,
    9.3677299909914948e+10,
    1.6819506294477527e+11,
    -1.6868283550187634e+12,
    6.9196677516269053e+12,
    -1.9436081834863168e+13,
    3.5195594052075742e+13,
    -2.8474916694776353e+12
  };
  double rc{0.0}, order{0.0};

  int ier = sixTerm(coeffs, scale, rc, order);

  ASSERT_THAT(ier, Eq(1));
  //EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  //EXPECT_THAT(order, DOUBLE_NEAR(3.1415E0));
  //EXPECT_THAT(coeffs.size(), Eq(30));
}

