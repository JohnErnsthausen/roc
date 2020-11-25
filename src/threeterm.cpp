#include <cfloat> // safemin
#include <cmath> // fabs
#include <string>
#include <vector>

extern "C"
{
#include "qrfactorization.h"
}
#include "data.hpp"
#include "exceptions.hpp"
#include "threeterm.hpp"

using namespace std;

// The three-term-test of Chang and Corliss
//
// The method returns |hrc - hrc_check|
//
// Whenever |hrc - hrc_check| > TOL the coefficients do not resemble a pole.
// 
// The calling subroutine should maintain the invarient that coeff.size()
// sufficiently large, say 10.
double threeterm(const vector<double> &coeff, const double &scale, double &rc,
              double &order)
{
  int nUse = THREETERM_NUSE;
  int m{nUse}, n{2}, k{0}, ier{0};
  string message;
  vector<int> ipiv(n, 0);
  vector<double> W(m*n, 0.0);
  vector<double> tau(n, 0.0);
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
  k = coeff.size() - nUse;
  for(int i{1}; i<=nUse; i++)
  {
    b(i) = k * coeff(k+1);
    W(i,1) = (k - 1) * coeff(k);
    W(i,2) = coeff(k);
    k++;
  }

  // Factor
  qrf(m, n, &W[0], m, &ipiv[0], &tau[0], &wrk[0], safmin, &ier);
  if ( ier != 0 )
  {
    message = "QRFactorization error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
  }
  // Solve
  qrs(m, n, &W[0], m, &tau[0], &b[0], &x[0], &ier);
  if ( ier != 0 )
  {
    message = "QRSolver error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
  }
  // Multiply by permutation
  for(int i = 1; i <=n; i++) { b(i) = x(ipiv(i)); }

  // Interpret the variables found from Least Squares Optimization Problem
  double hOverRc{0.0}, singularityOrder{0.0};

  hOverRc = b(1);
  rc = scale / hOverRc;
  if( isnan(rc) )
  {
    message = "Unconstrained optimization lead to NaN for Radius of Convergence\n";
    throw sayMessage(message);
  }

  singularityOrder = b(2) / b(1);
  if( isnan(singularityOrder) )
  {
    message = "Unconstrained optimization lead to NaN for Order of Singularity\n";
    throw sayMessage(message);
  }
  order = singularityOrder;

  // Evaluate previous W equation at the solution of the least squares optimization problem.
  // Return the backward error, the absolute value of this evaluation.
  nUse++;
  k = coeff.size() - nUse;
  double check = k * coeff(k+1) - ((k - 1) * coeff(k) * b(1) + coeff(k) * b(2));
  return fabs(check);
}
