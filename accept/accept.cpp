#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>
#include <cmath>
#include <numeric>
#include <vector>

#include "tadiff.h"
#include "data.hpp"
#include "matrix.hpp"
#include "vectorf.hpp"
#include "sixterm.hpp"

using namespace fadbad;
using namespace testing;

#define DIM 35
#define prec double

T<prec> sumPoles(const T<prec>& x, const int order)
{
  T<prec> z = 1. + 25.0 * x * x;
  z = 1.0/z;

  // x[0] = 1.0/(1.0+25.0*t*t) + 1.0/(1.0+25.0*t*t)^2;
  std::vector<T<prec>> powers(DIM, 1.0);
  for(int i{0}; i<DIM; i++)
  {
    for(int j{0}; j<=i; j++)
    {
      powers[j] *= z;
    }
  }
  reverse(powers.begin(), powers.end());
  return std::accumulate( powers.begin(), powers.begin()+order, T<prec> (0.0) );
}

TEST(AcceptanceTest,ComputePredictedPowerOf2)
{
  int orderActual{2};
  std::vector<double> coeffs(DIM);
  double a{1.0 / 5.0}, time{-0.3}, scale{1.0};

  T<prec> x, f;
  x = time; x[ 1 ] = scale; f = sumPoles(x,orderActual); f.eval(DIM);
  for (int i{0}; i < DIM; i++) { coeffs[i] = (double)f[i]; }
  
  double rc{0.0}, orderComputed{0.0};

  double err = sixterm(coeffs, scale, rc, orderComputed);

  EXPECT_THAT(err, DoubleNear(0.0, DBL_EPSILON));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), DBL_EPSILON));
  EXPECT_THAT(orderComputed, DoubleNear((double)orderActual, DBL_EPSILON));
  EXPECT_THAT(coeffs.size(), Eq(DIM));
}

TEST(AcceptanceTest,CorrectSumPowerTo3)
{
  int orderActual{3};
  std::vector<double> coeffs(DIM);
  double a{1.0 / 5.0}, time{-0.3}, scale{1.0};

  T<prec> x, f;
  x = time; x[ 1 ] = scale; f = sumPoles(x,orderActual); f.eval(DIM);
  for (int i{0}; i < DIM; i++) { coeffs[i] = (double)f[i]; }
  
  double rc{0.0}, orderComputed{0.0};

  double err = sixterm(coeffs, scale, rc, orderComputed);

  EXPECT_THAT(err, DoubleNear(0.0, DBL_EPSILON));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), DBL_EPSILON));
  EXPECT_THAT(orderComputed, DoubleNear((double)orderActual, DBL_EPSILON));
  EXPECT_THAT(coeffs.size(), Eq(DIM));
}

TEST(AcceptanceTest,CorrectSumPowerTo4)
{
  int orderActual{4};
  std::vector<double> coeffs(DIM);
  double a{1.0 / 5.0}, time{-0.3}, scale{1.0};

  T<prec> x, f;
  x = time; x[ 1 ] = scale; f = sumPoles(x,orderActual); f.eval(DIM);
  for (int i{0}; i < DIM; i++) { coeffs[i] = (double)f[i]; }
  
  double rc{0.0}, orderComputed{0.0};

  double err = sixterm(coeffs, scale, rc, orderComputed);

  EXPECT_THAT(err, DoubleNear(0.0, DBL_EPSILON));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), DBL_EPSILON));
  EXPECT_THAT(orderComputed, DoubleNear((double)orderActual, DBL_EPSILON));
  EXPECT_THAT(coeffs.size(), Eq(DIM));
}

TEST(AcceptanceTest,CorrectSumPowerTo7)
{
  int orderActual{7};
  std::vector<double> coeffs(DIM);
  double a{1.0 / 5.0}, time{-0.3}, scale{1.0};

  T<prec> x, f;
  x = time; x[ 1 ] = scale; f = sumPoles(x,orderActual); f.eval(DIM);
  for (int i{0}; i < DIM; i++) { coeffs[i] = (double)f[i]; }
  
  double rc{0.0}, orderComputed{0.0};

  double err = sixterm(coeffs, scale, rc, orderComputed);

  EXPECT_THAT(err, DoubleNear(0.0, DBL_EPSILON));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), DBL_EPSILON));
  EXPECT_THAT(orderComputed, DoubleNear((double)orderActual, DBL_EPSILON));
  EXPECT_THAT(coeffs.size(), Eq(DIM));
}

TEST(AcceptanceTest,Power8)
{
  double a{1.0 / 5.0}, time{-0.3}, scale{1.0}, rc{0.0}, order{0.0};
  
  std::vector<double> coeffs{
    8.0340237670171067e-05,
    2.9664087755140087e-03,
    5.6666013788665037e-02,
    7.4247982961090264e-01,
    7.4622665908384347e+00,
    6.0955696448353223e+01,
    4.1821679548767497e+02,
    2.4537928095396078e+03,
    1.2387851716468725e+04,
    5.3406958853580109e+04,
    1.9034041470099834e+05,
    5.0384790729501040e+05,
    5.0962523137455038e+05,
    -4.4310004262738079e+06,
    -3.8516544971975826e+07,
    -1.9483045193047231e+08,
    -7.3709802281675160e+08,
    -2.0698948047951820e+09,
    -3.1885834791594882e+09,
    7.5160430306028824e+09,
    8.8527590536210022e+10,
    4.4842564650093274e+11,
    1.6138506957170266e+12,
    4.1664092313235464e+12,
    5.1823350195233555e+12,
    -1.9381270051743824e+13,
    -1.7486465913941400e+14,
    -7.8991651698633412e+14,
    -2.5395415312674520e+15,
    -5.5405004901728730e+15,
    -2.8869957732739770e+15,
    4.5533255184943808e+16,
    2.8804807134809782e+17,
    1.1126085590967064e+18,
    3.0642224275153628e+18
  };

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, DBL_EPSILON));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), DBL_EPSILON));
  EXPECT_THAT(order, DoubleNear(8.0, DBL_EPSILON));
  EXPECT_THAT(coeffs.size(), Eq(DIM));
}

TEST(AcceptanceTest,SumTo7)
{
  double a{1.0 / 5.0}, time{-0.3}, scale{1.0}, rc{0.0}, order{0.0};
  
  std::vector<double> coeffs{
    4.4432839743447650e-01,
    2.9584401051385694e+00,
    1.4722984320930852e+01,
    6.4554014545202321e+01,
    2.6031279592621587e+02,
    9.7160336219037117e+02,
    3.3043576583642521e+03,
    9.7999172048031251e+03,
    2.2435898129030618e+04,
    1.8862613799904000e+04,
    -1.8006955738373377e+05,
    -1.4623031801268444e+06,
    -7.2399477899365267e+06,
    -2.7818980076471418e+07,
    -8.4583498851897508e+07,
    -1.8091207822478792e+08,
    -7.2474792231550753e+07,
    1.7676129529038100e+09,
    1.1585772847301514e+10,
    4.8483443100035431e+10,
    1.5215461712398621e+11,
    3.3460560468190338e+11,
    2.1484391780291943e+11,
    -2.5271964799268604e+12,
    -1.6876339061477297e+13,
    -6.8178334125683867e+13,
    -2.0121503091138594e+14,
    -3.9332291847646444e+14,
    -4.1373024550857125e+13,
    3.9439529510656000e+15,
    2.2203283141522488e+16,
    8.0744408732831824e+16,
    2.1115624119763402e+17,
    3.1784234369057779e+17,
    -4.3647429350057325e+17
  };

  double err = sixterm(coeffs, scale, rc, order);

  EXPECT_THAT(err, DoubleNear(0.0, DBL_EPSILON));
  EXPECT_THAT(rc, DoubleNear(sqrt(a * a + time * time), DBL_EPSILON));
  EXPECT_THAT(order, DoubleNear(7.0, DBL_EPSILON));
  EXPECT_THAT(coeffs.size(), Eq(DIM));
}

TEST(AcceptanceTest,At8)
{
  std::vector<double> coeffsAt8{
    8.0340237670171067e-05,
    2.9664087755140087e-03,
    5.6666013788665037e-02,
    7.4247982961090264e-01,
    7.4622665908384347e+00,
    6.0955696448353223e+01,
    4.1821679548767497e+02,
    2.4537928095396078e+03,
    1.2387851716468725e+04,
    5.3406958853580109e+04,
    1.9034041470099834e+05,
    5.0384790729501040e+05,
    5.0962523137455038e+05,
    -4.4310004262738079e+06,
    -3.8516544971975826e+07,
    -1.9483045193047231e+08,
    -7.3709802281675160e+08,
    -2.0698948047951820e+09,
    -3.1885834791594882e+09,
    7.5160430306028824e+09,
    8.8527590536210022e+10,
    4.4842564650093274e+11,
    1.6138506957170266e+12,
    4.1664092313235464e+12,
    5.1823350195233555e+12,
    -1.9381270051743824e+13,
    -1.7486465913941400e+14,
    -7.8991651698633412e+14,
    -2.5395415312674520e+15,
    -5.5405004901728730e+15,
    -2.8869957732739770e+15,
    4.5533255184943808e+16,
    2.8804807134809782e+17,
    1.1126085590967064e+18,
    3.0642224275153628e+18
  };

  std::vector<double> coeffsTo7{
    4.4432839743447650e-01,
    2.9584401051385694e+00,
    1.4722984320930852e+01,
    6.4554014545202321e+01,
    2.6031279592621587e+02,
    9.7160336219037117e+02,
    3.3043576583642521e+03,
    9.7999172048031251e+03,
    2.2435898129030618e+04,
    1.8862613799904000e+04,
    -1.8006955738373377e+05,
    -1.4623031801268444e+06,
    -7.2399477899365267e+06,
    -2.7818980076471418e+07,
    -8.4583498851897508e+07,
    -1.8091207822478792e+08,
    -7.2474792231550753e+07,
    1.7676129529038100e+09,
    1.1585772847301514e+10,
    4.8483443100035431e+10,
    1.5215461712398621e+11,
    3.3460560468190338e+11,
    2.1484391780291943e+11,
    -2.5271964799268604e+12,
    -1.6876339061477297e+13,
    -6.8178334125683867e+13,
    -2.0121503091138594e+14,
    -3.9332291847646444e+14,
    -4.1373024550857125e+13,
    3.9439529510656000e+15,
    2.2203283141522488e+16,
    8.0744408732831824e+16,
    2.1115624119763402e+17,
    3.1784234369057779e+17,
    -4.3647429350057325e+17
  };

  std::vector<double> coeffsTo8{
    4.4440873767214667e-01,
    2.9614065139140835e+00,
    1.4779650334719516e+01,
    6.5296494374813221e+01,
    2.6777506251705432e+02,
    1.0325590586387245e+03,
    3.7225744538519270e+03,
    1.2253710014342743e+04,
    3.4823749845499566e+04,
    7.2269572653485811e+04,
    1.0270857317250571e+04,
    -9.5845527283237898e+05,
    -6.7303225585703012e+06,
    -3.2249980502835833e+07,
    -1.2310004382466999e+08,
    -3.7574253016121304e+08,
    -8.0957281508707523e+08,
    -3.0228185211446023e+08,
    8.3971893670035362e+09,
    5.5999486125523071e+10,
    2.4068220764053662e+11,
    7.8303125112365405e+11,
    1.8286946134290198e+12,
    1.6392127518206421e+12,
    -1.1694004036890205e+13,
    -8.7559604146016312e+13,
    -3.7607968990134625e+14,
    -1.1832394348757695e+15,
    -2.5809145538934595e+15,
    -1.5965475339821690e+15,
    1.9316287378192748e+16,
    1.2627766392470070e+17,
    4.9920431249640806e+17,
    1.4304509024674401e+18,
    2.6277481326464077e+18
  };

  for(int i{0}; i<DIM; i++)
  {
    EXPECT_THAT(fabs((coeffsTo7[i]+coeffsAt8[i]) - coeffsTo8[i])/fabs(coeffsTo8[i]), DoubleNear(0.0, DBL_EPSILON));
  }
}
  
