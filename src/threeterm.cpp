#include <cmath>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "matrix.hpp"
#include "qrfactorization.hpp"
#include "vectorf.hpp"

void constructThreeTermSystem(const vectorf<double> &coeff, const int nUse,
                              matrix<double> &W, vectorf<double> &b)
{
  // TODO Check coeff.get_size, nUse, and W.get_rows are compatible?

  for (int i{1}, k{(int)coeff.get_size() - nUse}; i <= nUse; i++, k++)
  {
    W(i, 1) = (k - 1) * coeff(k);
    W(i, 2) = coeff(k);
    b(i) = k * coeff(k + 1);
  }
}

void testRCThree(double rc)
{
  if (std::isnan(rc))
  {
    std::string message =
        "Unconstrained optimization lead to NaN for Radius of Convergence\n";
    throw sayMessage(message);
  }
}

void testOrder(double order)
{
  if (std::isnan(order))
  {
    std::string message =
        "Unconstrained optimization lead to NaN for Order of Singularity\n";
    throw sayMessage(message);
  }
}

// The three-term-test of Chang and Corliss
double threeterm(const std::vector<double> &coeff, const double &scale,
                 double &rc, double &order)
{
  int nUse = THREETERM_NUSE;
  int m{nUse}, n{2};
  matrix<double> W(m, n);
  vectorf<double> beta(n);
  vectorf<double> b(m);
  vectorf<double> tc(coeff);

  constructThreeTermSystem(tc, nUse, W, b);

  // Solve W beta = b for beta
  qr(m, n, W, b, beta);

  // Interpret the variables found from Least Squares Optimization Problem
  double hOverRc = beta(1);
  rc = scale / hOverRc;
  testRCThree(rc);

  order = beta(2) / beta(1);
  testOrder(order);

  // Evaluate previous W equation at the solution of the least squares
  // optimization problem. Return the backward error, the absolute value of this
  // evaluation.
  //
  // TODO Use all equations?
  nUse++;
  int k = coeff.size() - nUse;
  double check = k * tc(k + 1) - ((k - 1) * tc(k) * beta(1) + tc(k) * beta(2));

  return fabs(check);
}
