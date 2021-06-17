#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <cfloat>
#include <cmath>

#include "sixterm.hpp"

#define epsilon DBL_EPSILON
using namespace testing;

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithScalingTenthForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t =-0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0 / 5.0};
  double time{-0.3};
  double scale{0.1};
  std::vector<double> coeffs{3.0769230769230771e-01,  1.4201183431952663e-01,
                        4.1875284478834776e-02,  8.4030671194986195e-03,
                        6.5716294139668780e-04,  -3.4308380547065328e-04,
                        -2.0889736724773902e-04, -7.0023107539675438e-05,
                        -1.6249329076177969e-05, -2.1132974551840289e-06,
                        2.7458033423644587e-07,  2.8929072773866954e-07,
                        1.1239723324581319e-07,  2.9622513210477662e-08,
                        5.0259881551579064e-09,  4.1031978497675339e-11,
                        -3.6767663724398882e-10, -1.7285321553550831e-10,
                        -5.1495588920697002e-11, -1.0470793691436441e-11,
                        -8.7147486368628090e-13, 4.0322650071682737e-13,
                        2.5314106676824962e-13,  8.5816915376359253e-14,
                        2.0135417345377379e-14,  2.6919683612234631e-15,
                        -3.0643132138743080e-16, -3.4850432996523446e-16,
                        -1.3727651218492121e-16, -3.6550364857253289e-17};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 5.12424e-14));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 1.88739e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 6.61694e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithScalingTenthForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0 / 5.0};
  double time{0.0};
  double scale{0.1};
  std::vector<double> coeffs{
      1.0000000000000000e+00,  0.0000000000000000e+00,  -2.5000000000000000e-01,
      0.0000000000000000e+00,  6.2500000000000000e-02,  0.0000000000000000e+00,
      -1.5625000000000000e-02, 0.0000000000000000e+00,  3.9062500000000000e-03,
      0.0000000000000000e+00,  -9.7656250000000000e-04, 0.0000000000000000e+00,
      2.4414062500000000e-04,  0.0000000000000000e+00,  -6.1035156250000000e-05,
      0.0000000000000000e+00,  1.5258789062500000e-05,  0.0000000000000000e+00,
      -3.8146972656250000e-06, 0.0000000000000000e+00,  9.5367431640625000e-07,
      0.0000000000000000e+00,  -2.3841857910156250e-07, 0.0000000000000000e+00,
      5.9604644775390625e-08,  0.0000000000000000e+00,  -1.4901161193847656e-08,
      0.0000000000000000e+00,  3.7252902984619141e-09,  0.0000000000000000e+00};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 2.60419e-15));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), epsilon));
  EXPECT_THAT(order, DoubleNear(1.0, 2.24266e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithScalingTenthForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0 / 5.0};
  double time{0.3};
  double scale{0.1};
  std::vector<double> coeffs{
      3.0769230769230771e-01,  -1.4201183431952663e-01, 4.1875284478834776e-02,
      -8.4030671194986195e-03, 6.5716294139668780e-04,  3.4308380547065328e-04,
      -2.0889736724773902e-04, 7.0023107539675438e-05,  -1.6249329076177969e-05,
      2.1132974551840289e-06,  2.7458033423644587e-07,  -2.8929072773866954e-07,
      1.1239723324581319e-07,  -2.9622513210477662e-08, 5.0259881551579064e-09,
      -4.1031978497675339e-11, -3.6767663724398882e-10, 1.7285321553550831e-10,
      -5.1495588920697002e-11, 1.0470793691436441e-11,  -8.7147486368628090e-13,
      -4.0322650071682737e-13, 2.5314106676824962e-13,  -8.5816915376359253e-14,
      2.0135417345377379e-14,  -2.6919683612234631e-15, -3.0643132138743080e-16,
      3.4850432996523446e-16,  -1.3727651218492121e-16, 3.6550364857253289e-17};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 5.12424e-14));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 1.88739e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 6.61694e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t =-0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0/5.0};
  double time{-0.3};
  double scale{1.0};
  std::vector<double> coeffs{3.0769230769230771e-01,  1.4201183431952664e+00,
                        4.1875284478834782e+00,  8.4030671194986173e+00,
                        6.5716294139668623e+00,  -3.4308380547065390e+01,
                        -2.0889736724773917e+02, -7.0023107539675470e+02,
                        -1.6249329076177976e+03, -2.1132974551840298e+03,
                        2.7458033423644597e+03,  2.8929072773866963e+04,
                        1.1239723324581322e+05,  2.9622513210477668e+05,
                        5.0259881551579072e+05,  4.1031978497674834e+04,
                        -3.6767663724398911e+06, -1.7285321553550843e+07,
                        -5.1495588920697041e+07, -1.0470793691436444e+08,
                        -8.7147486368627921e+07, 4.0322650071682835e+08,
                        2.5314106676824994e+09,  8.5816915376359320e+09,
                        2.0135417345377384e+10,  2.6919683612234612e+10,
                        -3.0643132138743221e+10, -3.4850432996523499e+11,
                        -1.3727651218492134e+12, -3.6550364857253306e+12};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 1.55432e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 1.34560e-13));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0/5.0};
  double time{0.0};
  double scale{1.0};
  std::vector<double> coeffs{
      1.0000000000000000e+00,  0.0000000000000000e+00,  -2.5000000000000000e+01,
      0.0000000000000000e+00,  6.2500000000000000e+02,  0.0000000000000000e+00,
      -1.5625000000000000e+04, 0.0000000000000000e+00,  3.9062500000000000e+05,
      0.0000000000000000e+00,  -9.7656250000000000e+06, 0.0000000000000000e+00,
      2.4414062500000000e+08,  0.0000000000000000e+00,  -6.1035156250000000e+09,
      0.0000000000000000e+00,  1.5258789062500000e+11,  0.0000000000000000e+00,
      -3.8146972656250000e+12, 0.0000000000000000e+00,  9.5367431640625000e+13,
      0.0000000000000000e+00,  -2.3841857910156250e+15, 0.0000000000000000e+00,
      5.9604644775390624e+16,  0.0000000000000000e+00,  -1.4901161193847657e+18,
      0.0000000000000000e+00,  3.7252902984619139e+19,  0.0000000000000000e+00};
  double rc{0.0}, order{0.0};

 double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), epsilon));
  EXPECT_THAT(order, DoubleNear(1.0, 4.44090e-16));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaOne)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  double a{1.0/5.0};
  double time{0.3};
  double scale{1.0};
  std::vector<double> coeffs{
      3.0769230769230771e-01,  -1.4201183431952664e+00, 4.1875284478834782e+00,
      -8.4030671194986173e+00, 6.5716294139668623e+00,  3.4308380547065390e+01,
      -2.0889736724773917e+02, 7.0023107539675470e+02,  -1.6249329076177976e+03,
      2.1132974551840298e+03,  2.7458033423644597e+03,  -2.8929072773866963e+04,
      1.1239723324581322e+05,  -2.9622513210477668e+05, 5.0259881551579072e+05,
      -4.1031978497674834e+04, -3.6767663724398911e+06, 1.7285321553550843e+07,
      -5.1495588920697041e+07, 1.0470793691436444e+08,  -8.7147486368627921e+07,
      -4.0322650071682835e+08, 2.5314106676824994e+09,  -8.5816915376359320e+09,
      2.0135417345377384e+10,  -2.6919683612234612e+10, -3.0643132138743221e+10,
      3.4850432996523499e+11,  -1.3727651218492134e+12, 3.6550364857253306e+12};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 1.55432e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 1.34560e-13));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithScalingHundredthForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t =-0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  double a{1.0 / 5.0};
  double time{-0.3};
  double scale{0.01};
  std::vector<double> coeffs{
      2.4655761562767709e-02,  3.5748957668969889e-03,  2.8208219200667809e-04,
      1.5636812703033154e-05,  6.5875429164410868e-07,  2.1094322169896426e-08,
      4.5261475074311576e-10,  1.1260909441922327e-12,  -4.6867806776404976e-13,
      -2.8056792096691283e-14, -1.0573048948187840e-15, -2.8313453349783412e-17,
      -4.3638374553989396e-19, 5.4965174709350763e-21,  7.3086407830397857e-22,
      3.3112669205272192e-23,  1.0201294809137287e-24,  2.1125410354692833e-26,
      1.1958555283498346e-28,  -1.3772058760031110e-29, -8.1538224313874953e-31,
      -2.8716158576770647e-32, -7.0504814826739304e-34, -9.3677299909588449e-36,
      1.6819506294660599e-37,  1.6868283550256514e-38,  6.9196677516458404e-40,
      1.9436081834894383e-41,  3.5195594052044087e-43,  2.8474916689612392e-46};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 3.38904e-13));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 3.02759e-13));
  EXPECT_THAT(order, DoubleNear(3.1415, 1.62773e-11));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithScalingThousandthForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  double a{1.0 / 5.0};
  double time{0.0};
  double scale{0.001};
  std::vector<double> coeffs{
      1.0000000000000000e+00,  0.0000000000000000e+00,  -7.8537500000000005e-05,
      0.0000000000000000e+00,  4.0657882031250008e-09,  0.0000000000000000e+00,
      -1.7420208371972664e-13, 0.0000000000000000e+00,  6.6866381072793826e-18,
      0.0000000000000000e+00,  -2.3876313021567852e-22, 0.0000000000000000e+00,
      8.0995417693789439e-27,  0.0000000000000000e+00,  -2.6443557530277716e-31,
      0.0000000000000000e+00,  8.3805418341659806e-36,  0.0000000000000000e+00,
      -2.5936613012600057e-40, 0.0000000000000000e+00,  7.8727346723120797e-45,
      0.0000000000000000e+00,  -2.3513532430952043e-49, 0.0000000000000000e+00,
      6.9274295598397387e-54,  0.0000000000000000e+00,  -2.0171475900060203e-58,
      0.0000000000000000e+00,  5.8142478257289360e-63,  0.0000000000000000e+00};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 1.11858e-14));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 1.58208e-15));
  EXPECT_THAT(order, DoubleNear(3.1415, 1.92292e-13));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithScalingThousandthForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  double a{1.0 / 5.0};
  double time{0.3};
  double scale{0.001};
  std::vector<double> coeffs{2.4655761562767709e-02,  -3.5748957668969888e-04,
                        2.8208219200667818e-06,  -1.5636812703033177e-08,
                        6.5875429164411143e-11,  -2.1094322169896626e-13,
                        4.5261475074312502e-16,  -1.1260909441924017e-19,
                        -4.6867806776406748e-21, 2.8056792096693676e-23,
                        -1.0573048948189600e-25, 2.8313453349792976e-28,
                        -4.3638374554030017e-31, -5.4965174709220135e-34,
                        7.3086407830372582e-36,  -3.3112669205275378e-38,
                        1.0201294809143237e-40,  -2.1125410354727341e-43,
                        1.1958555283637769e-46,  1.3772058759989342e-48,
                        -8.1538224313796743e-51, 2.8716158576775504e-53,
                        -7.0504814826864239e-56, 9.3677299910270096e-59,
                        1.6819506294408271e-61,  -1.6868283550188946e-63,
                        6.9196677516356937e-66,  -1.9436081834912157e-68,
                        3.5195594052256001e-71,  -2.8474916699531912e-75};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 2.2289e-13));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 6.43930e-14));
  EXPECT_THAT(order, DoubleNear(3.1415, 1.64136e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativePT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  double a{1.0/5.0};
  double time{-0.3};
  double scale{1.0};
  std::vector<double> coeffs{
      2.4655761562767709e-02,  3.5748957668969888e-01,  2.8208219200667819e+00,
      1.5636812703033181e+01,  6.5875429164411216e+01,  2.1094322169896708e+02,
      4.5261475074313199e+02,  1.1260909441928604e+02,  -4.6867806776404268e+03,
      -2.8056792096692519e+04, -1.0573048948189148e+05, -2.8313453349791735e+05,
      -4.3638374554030015e+05, 5.4965174709191988e+05,  7.3086407830349784e+06,
      3.3112669205263112e+07,  1.0201294809138143e+08,  2.1125410354710954e+08,
      1.1958555283602265e+08,  -1.3772058759989672e+09, -8.1538224313751421e+09,
      -2.8716158576746910e+10, -7.0504814826746704e+10, -9.3677299909914948e+10,
      1.6819506294477527e+11,  1.6868283550187634e+12,  6.9196677516269053e+12,
      1.9436081834863168e+13,  3.5195594052075742e+13,  2.8474916694776353e+12};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 2.38855e-15));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 5.38459e-14));
  EXPECT_THAT(order, DoubleNear(3.1415, 3.79919e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtZeroWithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  double a{1.0/5.0};
  double time{0.0};
  double scale{1.0};
  std::vector<double> coeffs{
      1.0000000000000000e+00,  0.0000000000000000e+00,  -7.8537500000000009e+01,
      0.0000000000000000e+00,  4.0657882031250006e+03,  0.0000000000000000e+00,
      -1.7420208371972659e+05, 0.0000000000000000e+00,  6.6866381072793789e+06,
      0.0000000000000000e+00,  -2.3876313021567836e+08, 0.0000000000000000e+00,
      8.0995417693789396e+09,  0.0000000000000000e+00,  -2.6443557530277716e+11,
      0.0000000000000000e+00,  8.3805418341659883e+12,  0.0000000000000000e+00,
      -2.5936613012600116e+14, 0.0000000000000000e+00,  7.8727346723121130e+15,
      0.0000000000000000e+00,  -2.3513532430952214e+17, 0.0000000000000000e+00,
      6.9274295598398177e+18,  0.0000000000000000e+00,  -2.0171475900060536e+20,
      0.0000000000000000e+00,  5.8142478257290687e+21,  0.0000000000000000e+00};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 7.21646e-16));
  EXPECT_THAT(order, DoubleNear(3.1415, 7.90489e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtPT3WithoutScalingForComplexConjugatePolesAtPMOneFifthAlphaThree)
{
  // t = 0.3;
  // x[0] = 1.0/(1.0+25.0*t*t)^{3.1415};

  double a{1.0/5.0};
  double time{0.3};
  double scale{1.0};
  std::vector<double> coeffs{2.4655761562767709e-02,  -3.5748957668969888e-01,
                        2.8208219200667819e+00,  -1.5636812703033181e+01,
                        6.5875429164411216e+01,  -2.1094322169896708e+02,
                        4.5261475074313199e+02,  -1.1260909441928604e+02,
                        -4.6867806776404268e+03, 2.8056792096692519e+04,
                        -1.0573048948189148e+05, 2.8313453349791735e+05,
                        -4.3638374554030015e+05, -5.4965174709191988e+05,
                        7.3086407830349784e+06,  -3.3112669205263112e+07,
                        1.0201294809138143e+08,  -2.1125410354710954e+08,
                        1.1958555283602265e+08,  1.3772058759989672e+09,
                        -8.1538224313751421e+09, 2.8716158576746910e+10,
                        -7.0504814826746704e+10, 9.3677299909914948e+10,
                        1.6819506294477527e+11,  -1.6868283550187634e+12,
                        6.9196677516269053e+12,  -1.9436081834863168e+13,
                        3.5195594052075742e+13,  -2.8474916694776353e+12};
  double rc{0.0}, order{0.0};

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 2.38855e-15));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), 5.38459e-14));
  EXPECT_THAT(order, DoubleNear(3.1415, 3.79919e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

// Test exception throws (Seems hard)
// QRFactorization and QRSolve tested in the developement of these methods
// SQRT tested above (by accident)
// The radius of convergence is infinity, which is highly unlikely (NOT TESTED)
// Unconstrained optimization lead to infinite CosTheta which is not in [-1, 1]
// (NOT TESTED) Unconstrained optimization lead to CosTheta [" +
// to_string(cosTheta) + "] not in [-1, 1] (by accident) Unconstrained
// optimization lead to NaN for Order of Singularity (NOT TESTED)
//
// Test real pole, which should fail.
TEST(
    SixTermAnalysisOf,
    TaylorSeriesAtNegativeOneNearARealPoleAtOneWithScalingTenthAlphaOneExpectedToFail)
{
  // t =-1.0;
  // x[0] = 1.0/(1.0-t);

  // double a{1.0};
  // double time{-1.0};
  double scale{0.1};
  std::vector<double> coeffs{0.5,
                        0.025,
                        0.00125,
                        6.25e-05,
                        3.125e-06,
                        1.5625e-07,
                        7.8125e-09,
                        3.90625e-10,
                        1.953125e-11,
                        9.765625e-13,
                        4.8828125e-14,
                        2.44140625e-15,
                        1.220703125e-16,
                        6.103515625e-18,
                        3.0517578125e-19,
                        1.52587890625e-20,
                        7.62939453125001e-22,
                        3.814697265625e-23,
                        1.9073486328125e-24,
                        9.53674316406251e-26,
                        4.76837158203126e-27,
                        2.38418579101563e-28,
                        1.19209289550781e-29,
                        5.96046447753907e-31,
                        2.98023223876954e-32,
                        1.49011611938477e-33,
                        7.45058059692384e-35,
                        3.72529029846192e-36,
                        1.86264514923096e-37,
                        9.3132257461548e-39};
  double rc{0.0}, order{0.0};

  // Unconstrained optimization lead to Sqrt of negative number: -0.360944
  ASSERT_THROW(sixterm(coeffs, scale, rc, order), std::exception);
}

TEST(SixTermAnalysisOf, TestBeta4)
{
  double beta4 = -1.0;
  EXPECT_THROW(testBeta4(beta4), std::exception);
}

TEST(SixTermAnalysisOf, TestRCSix)
{
  double rc = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(testRCSix(rc), std::exception);
}

TEST(SixTermAnalysisOf, TestCosThetaNAN)
{
  double cosTheta = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(testCosTheta(cosTheta), std::exception);
}

TEST(SixTermAnalysisOf, TestCosThetaLessThanNegativeOne)
{
  double cosTheta = -2.0;
  EXPECT_THROW(testCosTheta(cosTheta), std::exception);
}

TEST(SixTermAnalysisOf, TestCosThetaGreaterThanPositiveOne)
{
  double cosTheta = 2.0;
  EXPECT_THROW(testCosTheta(cosTheta), std::exception);
}

TEST(SixTermAnalysisOf, TestSingularityOrderBothNAN)
{
  double singularityOrder1 = std::numeric_limits<double>::quiet_NaN();
  double singularityOrder2 = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(testSingularityOrder(singularityOrder1, singularityOrder2),
               std::exception);
}

TEST(SixTermAnalysisOf,
     TestSingularityOrderSingularityOrderOneIsNANSingularityOrderTwoIsOne)
{
  double singularityOrder1 = std::numeric_limits<double>::quiet_NaN();
  double singularityOrder2 = 1.0;

  EXPECT_THAT(testSingularityOrder(singularityOrder1, singularityOrder2),
              DoubleNear(1.0, epsilon));
}

TEST(SixTermAnalysisOf,
     TestSingularityOrderSingularityOrderOneIsOneSingularityOrderTwoIsNAN)
{
  double singularityOrder1 = 1.0;
  double singularityOrder2 = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THAT(testSingularityOrder(singularityOrder1, singularityOrder2),
              DoubleNear(1.0, epsilon));
}

TEST(
    SixTermAnalysisOf,
    TestSingularityOrderSingularityOrderOneIsOneSingularityOrderTwoIsThreeOrderIsAverage)
{
  double singularityOrder1 = 1.0;
  double singularityOrder2 = 3.0;

  EXPECT_THAT(testSingularityOrder(singularityOrder1, singularityOrder2),
              DoubleNear(2.0, epsilon));
}

