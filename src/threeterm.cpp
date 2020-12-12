#include <cmath>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "qrfactorization.hpp"
#include "matrix.hpp"
#include "qrfactorization.hpp"
#include "vectorf.hpp"

using namespace std;

// The three-term-test of Chang and Corliss
double threeterm(const vector<double> &coeff, const double &scale, double &rc,
                 double &order)
{
  int nUse = THREETERM_NUSE;
  int m{nUse}, n{2};
  matrix<double> W(m, n);
  vectorf<double> beta(n);
  vectorf<double> b(m);
  vectorf<double> tc(coeff);

  for (int i{1}, k{(int)coeff.size() - nUse}; i <= nUse; i++, k++)
  {
    W(i, 1) = (k - 1) * tc(k);
    W(i, 2) = tc(k);
    b(i) = k * tc(k+1);
  }

  // Solve W beta = b for beta
  qr(m, n, W, b, beta);

  // Interpret the variables found from Least Squares Optimization Problem
  double hOverRc = beta(1);
  rc = scale / hOverRc;
  if (isnan(rc))
  {
    std::string message =
        "Unconstrained optimization lead to NaN for Radius of Convergence\n";
    throw sayMessage(message);
  }

  order = beta(2) / beta(1);
  if (isnan(order))
  {
    std::string message =
        "Unconstrained optimization lead to NaN for Order of Singularity\n";
    throw sayMessage(message);
  }

  // Evaluate previous W equation at the solution of the least squares
  // optimization problem. Return the backward error, the absolute value of this
  // evaluation.
  nUse++;
  int k = coeff.size() - nUse;
  double check =
      k * tc(k+1) - ((k - 1) * tc(k) * beta(1) + tc(k) * beta(2));

  return fabs(check);
}
