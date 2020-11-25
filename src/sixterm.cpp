#include <cfloat>  // safemin
#include <cmath>   // fabs
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

// The six-term-test of Chang and Corliss
//
// The method returns the absolute value of the penultimate equation to nUse
// parameter, which must be at least 4, evaluated at the solution to the least
// squares optimal solution.
//
// The computation is successful whenever the value return is acceptably small,
// otherwise the algorithm is said to detect that the coefficients do not
// resemble pair of complex conjugate poles.
//
// The calling subroutine should maintain the invarient that coeff.size()
// must be sufficiently large, say 10.
double sixterm(const vector<double> &coeff, const double &scale, double &rc,
               double &order)
{
  int nUse = SIXTERM_NUSE;
  int m{nUse}, n{4}, k{0}, ier{0};
  string message;
  vector<int> ipiv(n, 0);
  vector<double> W(m * n, 0.0);
  vector<double> tau(n, 0.0);
  vector<double> wrk(m, 0.0);
  vector<double> x(m, 0.0);
  vector<double> b(m, 0.0);
  double safmin{DBL_MIN};

#define map1(i) (i) - 1
#define b(i) b[ map1(i) ]
#define x(i) x[ map1(i) ]
#define ipiv(i) ipiv[ map1(i) ]
#define coeff(i) coeff[ map1(i) ]
#define map(i, j) ((j)-1) * m + ((i)-1)
#define W(i, j) W[ map(i, j) ]
  k = coeff.size() - nUse;
  for (int i{1}; i <= nUse; i++)
  {
    b(i) = k * coeff(k + 1);
    W(i, 1) = 2.0 * coeff(k);
    W(i, 2) = 2.0 * (k - 1) * coeff(k);
    W(i, 3) = -2.0 * coeff(k - 1);
    W(i, 4) = -(k - 2) * coeff(k - 1);
    k++;
  }

  // Factor
  qrf(m, n, &W[ 0 ], m, &ipiv[ 0 ], &tau[ 0 ], &wrk[ 0 ], safmin, &ier);
  if (ier != 0)
  {
    message = "QRFactorization error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
  }
  // Solve
  qrs(m, n, &W[ 0 ], m, &tau[ 0 ], &b[ 0 ], &x[ 0 ], &ier);
  if (ier != 0)
  {
    message = "QRSolver error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
  }
  // Multiply by permutation
  for (int i = 1; i <= n; i++)
  {
    b(i) = x(ipiv(i));
  }

  // Interpret the variables found from Least Squares Optimization Problem
  double hOverRc{0.0}, cosTheta{0.0}, singularityOrder1{0.0},
      singularityOrder2{0.0};

  // Evaluate h/Rc
  if (b(4) < 0)
  {
    message = "Unconstrained optimization lead to Sqrt of negative number: " +
              to_string(b(4)) + " \n";
    throw sayMessage(message);
  }
  hOverRc = sqrt(b(4));

  // Evaluate Rc
  rc = scale / hOverRc;
  if (isnan(rc))
  {
    message =
        "The radius of convergence is infinity, which is highly unlikely\n";
    throw sayMessage(message);
  }

  // Evaluate Cos(Theta)
  cosTheta = b(2) / hOverRc;
  if (isnan(cosTheta))
  {
    message =
        "Unconstrained optimization lead to infinite CosTheta which is not in "
        "[-1, 1]\n";
    throw sayMessage(message);
  }
  if ((cosTheta < -1.0) || (cosTheta > 1.0))
  {
    message = "Unconstrained optimization lead to CosTheta [" +
              to_string(cosTheta) + "] not in [-1, 1]\n";
    throw sayMessage(message);
  }

  // Evaluate Order of the Singularity
  singularityOrder1 = b(1) / b(2);
  singularityOrder2 = b(3) / b(4);
  if (isnan(singularityOrder1) && isnan(singularityOrder2))
  {
    message =
        "Unconstrained optimization lead to NaN for Order of Singularity\n";
    throw sayMessage(message);
  }

  if (isnan(singularityOrder1) && !isnan(singularityOrder2))
    order = singularityOrder2;

  if (!isnan(singularityOrder1) && isnan(singularityOrder2))
    order = singularityOrder1;

  if (!isnan(singularityOrder1) && !isnan(singularityOrder2))
    order = (singularityOrder1 + singularityOrder2) / 2.0;

  // Compare order
  // printf( "Abs of difference between two computations for order:
  // [%22.16f]\n",
  //     fabs(singularityOrder1-singularityOrder2));
  // TODO Add this error to checked error on return?

  // Evaluate previous W equation at the solution of the least squares
  // optimization problem. Return the backward error, the absolute value of this
  // evaluation.
  nUse++;
  k = coeff.size() - nUse;
  double check =
      k * coeff(k + 1) -
      ((2.0 * coeff(k)) * b(1) + (2.0 * (k - 1) * coeff(k)) * b(2) +
       (-2.0 * coeff(k - 1)) * b(3) + (-(k - 2) * coeff(k - 1)) * b(4));

  return fabs(check);
}
