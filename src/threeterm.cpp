#include <cmath>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "qrfactorization.hpp"
#include "threeterm.hpp"

using namespace std;

// The three-term-test of Chang and Corliss
double threeterm(const vector<double> &coeff, const double &scale, double &rc,
                 double &order)
{
  int nUse = THREETERM_NUSE;
  int m{nUse}, n{2}, k{0};
  string message;
  vector<double> W(m * n, 0.0);
  vector<double> beta(n, 0.0);
  vector<double> b(m, 0.0);

#define map1(i) (i) - 1
#define b(i) b[ map1(i) ]
#define beta(i) beta[ map1(i) ]
#define coeff(i) coeff[ map1(i) ]
#define map(i, j) ((j)-1) * m + ((i)-1)
#define W(i, j) W[ map(i, j) ]
  k = coeff.size() - nUse;
  for (int i{1}; i <= nUse; i++)
  {
    b(i) = k * coeff(k + 1);
    W(i, 1) = (k - 1) * coeff(k);
    W(i, 2) = coeff(k);
    k++;
  }

  // Solve W beta = b for beta
  qr(m, n, W, b, beta);

  // Interpret the variables found from Least Squares Optimization Problem
  double hOverRc = beta(1);
  rc = scale / hOverRc;
  if (isnan(rc))
  {
    message =
        "Unconstrained optimization lead to NaN for Radius of Convergence\n";
    throw sayMessage(message);
  }

  order = beta(2) / beta(1);
  if (isnan(order))
  {
    message =
        "Unconstrained optimization lead to NaN for Order of Singularity\n";
    throw sayMessage(message);
  }

  // Evaluate previous W equation at the solution of the least squares
  // optimization problem. Return the backward error, the absolute value of this
  // evaluation.
  nUse++;
  k = coeff.size() - nUse;
  double check =
      k * coeff(k + 1) - ((k - 1) * coeff(k) * beta(1) + coeff(k) * beta(2));

  return fabs(check);
}
