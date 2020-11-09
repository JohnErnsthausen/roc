#include <cmath>
#include <iostream>
#include <vector>
#include "threeterm.hpp"

using namespace std;

// The three-term-test of Chang and Corliss
//
// The vector coeff is expected to have at least length 10
//
// The method returns ier=0 if computation was successful. Otherwise
// the method returns ier=1 whenever the algorithm detects that
// the coefficients do not resemble a pole.
int threeTerm(const vector<double> &coeff, const double &scale, double &rc, double &order)
{
  const double zero{0.0};

  int k{0}, n{0}, n_check{0}, ier{0};
  double hrc{0.0}, hrc_check{0.0};
  double vL{0.0}, vM{0.0}, vN{0.0};
  double wL{0.0}, wM{0.0}, wN{0.0};

  // Check for sufficient data to perform three term analysis
  k = (int) coeff.size();
  if( k < MINTERMS ) throw morecoefficients();
  
  // Extract coefficients and check for divide by zero
  n_check = k-5;
  wL = coeff[k-6];
  wM = coeff[k-5];
  wN = coeff[k-4];
  if( wM == zero || wL == zero ) throw dividebyzero();
  
  n = k-2;
  vL = coeff[k-3];
  vM = coeff[k-2];
  vN = coeff[k-1];
  if( vM == zero || vL == zero ) throw dividebyzero();

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

