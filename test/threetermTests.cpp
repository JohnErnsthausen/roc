#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <cfloat>
#include <cmath>
#include <limits>

#include "threeterm.hpp"

#define epsilon DBL_EPSILON
using namespace testing;

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtNegativeOneRealPoleAtOneWithScalingTenthAlphaOne)
{
  // t =-1.0;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{-1.0};
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

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 6.26323e-14));
  EXPECT_THAT(rc, DoubleNear(a - time, 8.57093e-14));
  EXPECT_THAT(order, DoubleNear(1.0, 1.05228e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtZeroOneRealPoleAtOneWithScalingTenthAlphaOne)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{0.0};
  double scale{0.1};
  std::vector<double> coeffs{1,     0.1,   0.01,  0.001, 0.0001, 1e-05, 1e-06, 1e-07,
                        1e-08, 1e-09, 1e-10, 1e-11, 1e-12,  1e-13, 1e-14, 1e-15,
                        1e-16, 1e-17, 1e-18, 1e-19, 1e-20,  1e-21, 1e-22, 1e-23,
                        1e-24, 1e-25, 1e-26, 1e-27, 1e-28,  1e-29};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 1.03287e-15));
  EXPECT_THAT(rc, DoubleNear(a - time, 4.44090e-16));
  EXPECT_THAT(order, DoubleNear(1.0, 1.35448e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtPTNineOneRealPoleAtOneWithScalingTenthAlphaOne)
{
  // t = 0.9;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{0.9};
  double scale{0.1};
  std::vector<double> coeffs{10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001,
                        10.0000000000001};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 2.2994e-15));
  EXPECT_THAT(rc, DoubleNear(a - time, epsilon));
  EXPECT_THAT(order, DoubleNear(1.0, 4.78507e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtOnePT0001RealPoleAtOneWithScalingTenthAlphaOne)
{
  // t = 1.0001;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{1.0001};
  double scale{0.1};
  std::vector<double> coeffs{
      -10000.0000000011,     10000000.0000022,      -10000000000.0033,
      10000000000004.4,      -1.00000000000055e+16, 1.00000000000066e+19,
      -1.00000000000077e+22, 1.00000000000088e+25,  -1.00000000000099e+28,
      1.0000000000011e+31,   -1.00000000000121e+34, 1.00000000000132e+37,
      -1.00000000000143e+40, 1.00000000000154e+43,  -1.00000000000165e+46,
      1.00000000000176e+49,  -1.00000000000187e+52, 1.00000000000198e+55,
      -1.00000000000209e+58, 1.0000000000022e+61,   -1.00000000000231e+64,
      1.00000000000242e+67,  -1.00000000000253e+70, 1.00000000000264e+73,
      -1.00000000000275e+76, 1.00000000000286e+79,  -1.00000000000297e+82,
      1.00000000000308e+85,  -1.00000000000319e+88, 1.00000000000331e+91};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(a - time, epsilon));
  EXPECT_THAT(order, DoubleNear(1.0, 8.03569e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtOnePT9RealPoleAtOneWithScalingTenthAlphaOne)
{
  // t = 1.9;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{1.9};
  double scale{0.1};
  std::vector<double> coeffs{
      -1.11111111111111,     0.123456790123457,     -0.0137174211248285,
      0.00152415790275873,   -0.000169350878084303, 1.88167642315892e-05,
      -2.09075158128769e-06, 2.32305731254188e-07,  -2.5811747917132e-08,
      2.86797199079245e-09,  -3.18663554532494e-10, 3.54070616147216e-11,
      -3.93411795719128e-12, 4.37124217465698e-13,  -4.85693574961887e-14,
      5.3965952773543e-15,   -5.99621697483812e-16, 6.66246330537568e-17,
      -7.40273700597298e-18, 8.22526333996998e-19,  -9.13918148885554e-20,
      1.01546460987284e-20,  -1.12829401096982e-21, 1.25366001218869e-22,
      -1.39295556909854e-23, 1.54772841010949e-24,  -1.71969823345499e-25,
      1.91077581494999e-26,  -2.12308423883332e-27, 2.35898248759258e-28};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 1.13933e-13));
  EXPECT_THAT(rc, DoubleNear(a - time, 7.76047e-14));
  EXPECT_THAT(order, DoubleNear(1.0, 1.98886e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf, TaylorSeriesAtThreePoleAtOneWithScalingTenthAlphaOne)
{
  // t = 3.0;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{3.0};
  double scale{0.1};
  std::vector<double> coeffs{-0.5,
                        0.025,
                        -0.00125,
                        6.25e-05,
                        -3.125e-06,
                        1.5625e-07,
                        -7.8125e-09,
                        3.90625e-10,
                        -1.953125e-11,
                        9.765625e-13,
                        -4.8828125e-14,
                        2.44140625e-15,
                        -1.220703125e-16,
                        6.103515625e-18,
                        -3.0517578125e-19,
                        1.52587890625e-20,
                        -7.62939453125001e-22,
                        3.814697265625e-23,
                        -1.9073486328125e-24,
                        9.53674316406251e-26,
                        -4.76837158203126e-27,
                        2.38418579101563e-28,
                        -1.19209289550781e-29,
                        5.96046447753907e-31,
                        -2.98023223876954e-32,
                        1.49011611938477e-33,
                        -7.45058059692384e-35,
                        3.72529029846192e-36,
                        -1.86264514923096e-37,
                        9.3132257461548e-39};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 6.26323e-14));
  EXPECT_THAT(rc, DoubleNear(a - time, 8.57093e-14));
  EXPECT_THAT(order, DoubleNear(1.0, 1.05228e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtNegativeOneRealPoleAtOneWithScalingOneAlphaOne)
{
  // t =-1.0;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{-1.0};
  double scale{1.0};
  std::vector<double> coeffs{
      5.000000000000000e-01, 2.500000000000000e-01, 1.250000000000000e-01,
      6.250000000000000e-02, 3.125000000000000e-02, 1.562500000000000e-02,
      7.812500000000000e-03, 3.906250000000000e-03, 1.953125000000000e-03,
      9.765625000000000e-04, 4.882812500000000e-04, 2.441406250000000e-04,
      1.220703125000000e-04, 6.103515625000000e-05, 3.051757812500000e-05,
      1.525878906250000e-05, 7.629394531250000e-06, 3.814697265625000e-06,
      1.907348632812500e-06, 9.536743164062500e-07, 4.768371582031250e-07,
      2.384185791015625e-07, 1.192092895507813e-07, 5.960464477539063e-08,
      2.980232238769531e-08, 1.490116119384766e-08, 7.450580596923828e-09,
      3.725290298461914e-09, 1.862645149230957e-09, 9.313225746154785e-10};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 4.52325e-15));
  EXPECT_THAT(rc, DoubleNear(a - time, 7.10544e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 8.33778e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtZeroOneRealPoleAtOneWithScalingOneAlphaOne)
{
  // t = 0.0;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{0.0};
  double scale{1.0};
  std::vector<double> coeffs{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 6.15225e-16));
  EXPECT_THAT(rc, DoubleNear(a - time, 6.66135e-16));
  EXPECT_THAT(order, DoubleNear(1.0, 2.40919e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtPTNineOneRealPoleAtOneWithScalingOneAlphaOne)
{
  // t = 0.9;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{0.9};
  double scale{1.0};
  std::vector<double> coeffs{
      1.000000000000000e+01, 1.000000000000000e+02, 1.000000000000001e+03,
      1.000000000000001e+04, 1.000000000000001e+05, 1.000000000000001e+06,
      1.000000000000002e+07, 1.000000000000002e+08, 1.000000000000002e+09,
      1.000000000000002e+10, 1.000000000000002e+11, 1.000000000000003e+12,
      1.000000000000003e+13, 1.000000000000003e+14, 1.000000000000003e+15,
      1.000000000000004e+16, 1.000000000000004e+17, 1.000000000000004e+18,
      1.000000000000004e+19, 1.000000000000005e+20, 1.000000000000005e+21,
      1.000000000000005e+22, 1.000000000000005e+23, 1.000000000000005e+24,
      1.000000000000006e+25, 1.000000000000006e+26, 1.000000000000006e+27,
      1.000000000000006e+28, 1.000000000000007e+29, 1.000000000000007e+30};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(a - time, 2.85883e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 8.24009e-13));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtOnePT0001RealPoleAtOneWithScalingOneAlphaOne)
{
  // t = 1.0001;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{1.0001};
  double scale{1.0};
  std::vector<double> coeffs{
      -1.000000000000110e+04,  1.000000000000220e+08,   -1.000000000000330e+12,
      1.000000000000440e+16,   -1.000000000000551e+20,  1.000000000000661e+24,
      -1.000000000000771e+28,  1.000000000000881e+32,   -1.000000000000991e+36,
      1.000000000001101e+40,   -1.000000000001211e+44,  1.000000000001321e+48,
      -1.000000000001431e+52,  1.000000000001542e+56,   -1.000000000001652e+60,
      1.000000000001762e+64,   -1.000000000001872e+68,  1.000000000001982e+72,
      -1.000000000002092e+76,  1.000000000002202e+80,   -1.000000000002312e+84,
      1.000000000002422e+88,   -1.000000000002533e+92,  1.000000000002643e+96,
      -1.000000000002753e+100, 1.000000000002863e+104,  -1.000000000002973e+108,
      1.000000000003083e+112,  -1.000000000003193e+116, 1.000000000003303e+120};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(a - time, epsilon));
  EXPECT_THAT(order, DoubleNear(1.0, 1.5100e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtOnePT9RealPoleAtOneWithScalingOneAlphaOne)
{
  // t = 1.9;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{1.9};
  double scale{1.0};
  std::vector<double> coeffs{
      -1.111111111111111e+00, 1.234567901234568e+00,  -1.371742112482853e+00,
      1.524157902758726e+00,  -1.693508780843029e+00, 1.881676423158922e+00,
      -2.090751581287691e+00, 2.323057312541879e+00,  -2.581174791713199e+00,
      2.867971990792444e+00,  -3.186635545324938e+00, 3.540706161472154e+00,
      -3.934117957191282e+00, 4.371242174656981e+00,  -4.856935749618868e+00,
      5.396595277354298e+00,  -5.996216974838109e+00, 6.662463305375677e+00,
      -7.402737005972976e+00, 8.225263339969974e+00,  -9.139181488855527e+00,
      1.015464609872837e+01,  -1.128294010969818e+01, 1.253660012188687e+01,
      -1.392955569098542e+01, 1.547728410109491e+01,  -1.719698233454990e+01,
      1.910775814949989e+01,  -2.123084238833321e+01, 2.358982487592579e+01};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 6.94819e-16));
  EXPECT_THAT(rc, DoubleNear(a - time, 2.44250e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 7.72716e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf, TaylorSeriesAtThreePoleAtOneWithScalingOneAlphaOne)
{
  // t = 3.0;
  // x[0] = 1.0/(1.0-t);

  double a{1.0};
  double time{3.0};
  double scale{1.0};
  std::vector<double> coeffs{
      -5.000000000000000e-01, 2.500000000000000e-01,  -1.250000000000000e-01,
      6.250000000000000e-02,  -3.125000000000000e-02, 1.562500000000000e-02,
      -7.812500000000000e-03, 3.906250000000000e-03,  -1.953125000000000e-03,
      9.765625000000000e-04,  -4.882812500000000e-04, 2.441406250000000e-04,
      -1.220703125000000e-04, 6.103515625000000e-05,  -3.051757812500000e-05,
      1.525878906250000e-05,  -7.629394531250000e-06, 3.814697265625000e-06,
      -1.907348632812500e-06, 9.536743164062500e-07,  -4.768371582031250e-07,
      2.384185791015625e-07,  -1.192092895507813e-07, 5.960464477539063e-08,
      -2.980232238769531e-08, 1.490116119384766e-08,  -7.450580596923828e-09,
      3.725290298461914e-09,  -1.862645149230957e-09, 9.313225746154785e-10};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 4.52325e-15));
  EXPECT_THAT(rc, DoubleNear(a - time, 7.10544e-15));
  EXPECT_THAT(order, DoubleNear(1.0, 8.33778e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

// Requires TOL=1e-1
TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtNegativeOnePoleAtOneWithScalingTenthAlphaTwo)
{
  // t =-1.0;
  // x[0] = 1.0/((1.0-t)*(1.0-t));

  double a{1.0};
  double time{-1.0};
  double scale{0.1};
  std::vector<double> coeffs{0.25,
                        0.025,
                        0.001875,
                        0.000125,
                        7.8125e-06,
                        4.6875e-07,
                        2.734375e-08,
                        1.5625e-09,
                        8.78906250000001e-11,
                        4.8828125e-12,
                        2.685546875e-13,
                        1.46484375e-14,
                        7.9345703125e-16,
                        4.2724609375e-17,
                        2.288818359375e-18,
                        1.220703125e-19,
                        6.48498535156251e-21,
                        3.4332275390625e-22,
                        1.81198120117188e-23,
                        9.53674316406251e-25,
                        5.00679016113282e-26,
                        2.62260437011719e-27,
                        1.37090682983399e-28,
                        7.15255737304688e-30,
                        3.72529029846192e-31,
                        1.9371509552002e-32,
                        1.00582838058472e-33,
                        5.21540641784669e-35,
                        2.70083546638489e-36,
                        1.39698386192322e-37};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 1.29579e-13));
  EXPECT_THAT(rc, DoubleNear(a - time, 2.07613e-13));
  EXPECT_THAT(order, DoubleNear(2.0, 2.5696e-12));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf, TaylorSeriesAtZeroPoleAtOneWithScalingTenthAlphaTwo)
{
  // t = 0.0;
  // x[0] = 1.0/((1.0-t)*(1.0-t));

  double a{1.0};
  double time{0.0};
  double scale{0.1};
  std::vector<double> coeffs{1,       0.2,     0.03,
                        0.004,   0.0005,  6e-05,
                        7e-06,   8e-07,   9.00000000000001e-08,
                        1e-08,   1.1e-09, 1.2e-10,
                        1.3e-11, 1.4e-12, 1.5e-13,
                        1.6e-14, 1.7e-15, 1.8e-16,
                        1.9e-17, 2e-18,   2.1e-19,
                        2.2e-20, 2.3e-21, 2.4e-22,
                        2.5e-23, 2.6e-24, 2.7e-25,
                        2.8e-26, 2.9e-27, 3e-28};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 9.30465e-15));
  EXPECT_THAT(rc, DoubleNear(a - time, 7.32748e-15));
  EXPECT_THAT(order, DoubleNear(2.0, 1.83409e-13));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf, TaylorSeriesAtPT9PoleAtOneWithScalingTenthAlphaTwo)
{
  // t = 0.9;
  // x[0] = 1.0/((1.0-t)*(1.0-t));

  double a{1.0};
  double time{0.9};
  double scale{0.1};
  std::vector<double> coeffs{100,
                        200,
                        300,
                        400.000000000001,
                        500.000000000001,
                        600.000000000002,
                        700.000000000002,
                        800.000000000003,
                        900.000000000004,
                        1000,
                        1100.00000000001,
                        1200.00000000001,
                        1300.00000000001,
                        1400.00000000001,
                        1500.00000000001,
                        1600.00000000001,
                        1700.00000000001,
                        1800.00000000001,
                        1900.00000000002,
                        2000.00000000002,
                        2100.00000000002,
                        2200.00000000002,
                        2300.00000000003,
                        2400.00000000003,
                        2500.00000000003,
                        2600.00000000004,
                        2700.00000000004,
                        2800.00000000004,
                        2900.00000000005,
                        3000.00000000005};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 2.12876e-15));
  EXPECT_THAT(rc, DoubleNear(a - time, 2.77557e-16));
  EXPECT_THAT(order, DoubleNear(2.0, 7.59394e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf,
     TaylorSeriesAtNegativeOnePoleAtOneWithScalingOneAlphaTwo)
{
  // t =-1.0;
  // x[0] = 1.0/((1.0-t)*(1.0-t));

  double a{1.0};
  double time{-1.0};
  double scale{1.0};
  std::vector<double> coeffs{
      2.500000000000000e-01, 2.500000000000000e-01, 1.875000000000000e-01,
      1.250000000000000e-01, 7.812500000000000e-02, 4.687500000000000e-02,
      2.734375000000000e-02, 1.562500000000000e-02, 8.789062500000000e-03,
      4.882812500000000e-03, 2.685546875000000e-03, 1.464843750000000e-03,
      7.934570312500000e-04, 4.272460937500000e-04, 2.288818359375000e-04,
      1.220703125000000e-04, 6.484985351562500e-05, 3.433227539062500e-05,
      1.811981201171875e-05, 9.536743164062500e-06, 5.006790161132813e-06,
      2.622604370117188e-06, 1.370906829833984e-06, 7.152557373046875e-07,
      3.725290298461914e-07, 1.937150955200195e-07, 1.005828380584717e-07,
      5.215406417846680e-08, 2.700835466384888e-08, 1.396983861923218e-08};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 4.79343e-16));
  EXPECT_THAT(rc, DoubleNear(a - time, 1.33228e-15));
  EXPECT_THAT(order, DoubleNear(2.0, 1.28787e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf, TaylorSeriesAtZeroPoleAtOneWithScalingOneAlphaTwo)
{
  // t = 0.0;
  // x[0] = 1.0/((1.0-t)*(1.0-t));

  double a{1.0};
  double time{0.0};
  double scale{1.0};
  std::vector<double> coeffs{1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                        11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                        21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, 3.25648e-16));
  EXPECT_THAT(rc, DoubleNear(a - time, 6.66135e-16));
  EXPECT_THAT(order, DoubleNear(2.0, 2.26486e-14));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf, TaylorSeriesAtPT9PoleAtOneWithScalingOneAlphaTwo)
{
  // t = 0.9;
  // x[0] = 1.0/((1.0-t)*(1.0-t));

  double a{1.0};
  double time{0.9};
  double scale{1.0};
  std::vector<double> coeffs{
      1.000000000000001e+02, 2.000000000000002e+03, 3.000000000000004e+04,
      4.000000000000006e+05, 5.000000000000010e+06, 6.000000000000015e+07,
      7.000000000000021e+08, 8.000000000000028e+09, 9.000000000000035e+10,
      1.000000000000004e+12, 1.100000000000005e+13, 1.200000000000006e+14,
      1.300000000000008e+15, 1.400000000000009e+16, 1.500000000000010e+17,
      1.600000000000010e+18, 1.700000000000012e+19, 1.800000000000013e+20,
      1.900000000000014e+21, 2.000000000000016e+22, 2.100000000000018e+23,
      2.200000000000020e+24, 2.300000000000022e+25, 2.400000000000025e+26,
      2.500000000000028e+27, 2.600000000000031e+28, 2.700000000000035e+29,
      2.800000000000039e+30, 2.900000000000042e+31, 3.000000000000046e+32};
  double rc{0.0}, order{0.0};

  double err = threeterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, epsilon));
  EXPECT_THAT(rc, DoubleNear(a - time, 1.23513e-15));
  EXPECT_THAT(order, DoubleNear(2.0, 3.52164e-13));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(
    ThreeTermAnalysisOf,
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

  double err = threeterm(coeffs, scale, rc, order);

  // Model is wrong, as it should be. However there are no difficulties
  // executing this test case
  EXPECT_THAT(err, DoubleNear(0.0, 18.0446));
  EXPECT_THAT(rc, DoubleNear(a - time, 0.545748));
  EXPECT_THAT(order, DoubleNear(1.0, 26.5677));
  EXPECT_THAT(coeffs.size(), Eq(30));
}

TEST(ThreeTermAnalysisOf, TestRCThree)
{
  double rc = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(testRCThree(rc), std::exception);
}

TEST(ThreeTermAnalysisOf, TestOrder)
{
  double order = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(testOrder(order), std::exception);
}

// WARNING Rc and order can be accurate, but backward error can be large near
// singularity WARNING Error model wrong for complex-conjugate comparison
//
// Exceptions (Hard to think of a test)
// Top-line comparison
