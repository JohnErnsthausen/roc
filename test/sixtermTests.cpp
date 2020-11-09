#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cmath>
#include <vector>
#include <cfloat>

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
int sixTerm(const vector<double> &coeff, const double &scale, double &rc, double &order)
{
  const double zero{0.0};

  int k{0}, n{0}, n_check{0}, ier{0};
  double hrc{0.0}, hrc_check{0.0};
  double vL{0.0}, vM{0.0}, vN{0.0};
  double wL{0.0}, wM{0.0}, wN{0.0};

  // Check for sufficient data to perform three term analysis
  k = (int) coeff.size();
  //if( k < MINTERMS ) throw morecoefficients();
  
  // Extract coefficients and check for divide by zero
  n_check = k-5;
  wL = coeff[k-6];
  wM = coeff[k-5];
  wN = coeff[k-4];
  //if( wM == zero || wL == zero ) throw dividebyzero();
  
  n = k-2;
  vL = coeff[k-3];
  vM = coeff[k-2];
  vN = coeff[k-1];
  //if( vM == zero || vL == zero ) throw dividebyzero();

  // This is the three term analysis.
  // TODO Perhaps make a method for this formula?
  hrc = n*(vN/vM) - (n-1)*(vM/vL);
  hrc_check = n_check*(wN/wM) - (n_check-1)*(wM/wL);

  // Check for agreement between computations against TOL
  //
  // If no agreement, three term analysis is said to have failed
  // becuase the coefficients do not represent a pole.
  if( fabs(hrc - hrc_check) > TOL )
  {
    cout.precision(16);
    cout << scientific;
    cout << "|hrc - hrc_check| = " << fabs(hrc - hrc_check) << "\n";
    ier = 1;
    return ier;
  }

  rc  = 1.0/hrc;
  order = n*(vN/vM)*rc - (n - 1);
  rc *= scale;

  return ier;
}

TEST(SixTermAnalysisOf, TaylorSeriesAtNegativePT3ComplexConjugatePolesAtPMOneFifthIWithScalingTenthAlphaOne)
{
  // t =-0.3;
  // x[0] = 1.0/(1.0+25.0*t*t);

  //double a{1.0/5.0};
  double time{-0.3};
  double scale{0.1};
  vector<double> coeffs{
     0.307692307692308,
     0.142011834319527,
    0.0418752844788348,
   0.00840306711949862,
  0.000657162941396686,
 -0.000343083805470654,
 -0.000208897367247739,
 -7.00231075396755e-05,
  -1.6249329076178e-05,
 -2.11329745518403e-06,
  2.74580334236446e-07,
   2.8929072773867e-07,
  1.12397233245813e-07,
  2.96225132104777e-08,
  5.02598815515791e-09,
  4.1031978497680128e-11,   // 4.10319784976749e-11;
 -3.67676637243989e-10,
 -1.72853215535509e-10,
 -5.14955889206971e-11,
 -1.04707936914365e-11,
  -8.7147486368628e-13,
  4.03226500716829e-13,
   2.5314106676825e-13,
  8.58169153763594e-14,
  2.01354173453774e-14,
  2.69196836122347e-15,
 -3.06431321387433e-16,
 -3.48504329965235e-16,
 -1.37276512184922e-16,
 -3.65503648572534e-17 
  };
  double rc{0.0}, order{0.0};

  const double safmin = DBL_MIN;
  const int m{3}, n{3}, lda{3};
  int *ipiv{new int[3]()};
  double *x{new double[3]()};
  double *b{new double[3]()};
  double *tau{new double[3]()};
  double *wrk{new double[3]()};
  double *a{new double[9]()};
  
  a[0]= 2.0;
  a[1]=-1.0;
  a[2]= 0.0;
  a[3]=-1.0;
  a[4]= 2.0;
  a[5]=-1.0;
  a[6]= 0.0;
  a[7]=-1.0;
  a[8]= 2.0;
  b[0]=-1.0;
  b[1]=-1.0;
  b[2]=-1.0;

  int ier{3}, res{0};
  res = qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier);
  if (res != 0 || ier != 0)
  {
    printf(
        "Error finding the Least Squares Solution with res[%d]!=0 or "
        "ier[%d]!=0\n",
        res, ier);
  }
  
  res = qrs(m, n, a, lda, tau, b, x, &ier);

  if (res != 0 || ier != 0)
  {
    printf(
        "Error finding the Least Squares Solution with res[%d]!=0 or "
        "ier[%d]!=0\n",
        res, ier);
  }

  for(int i=0; i<3; i++)
    printf("b[%d] = [%22.16f]\n",i,x[ipiv[i]]);

  //  int ier = sixTerm(coeffs, scale, rc, order);

  //  EXPECT_THAT(ier, Eq(0));
  //  EXPECT_THAT(rc, DOUBLE_NEAR(sqrt(a*a+time*time)));
  //  EXPECT_THAT(order, DOUBLE_NEAR(1.0));
  //  EXPECT_THAT(coeffs.size(), Eq(30));
  delete [] ipiv;
  delete [] x;
  delete [] b;
  delete [] tau;
  delete [] wrk;
  delete [] a;
}
