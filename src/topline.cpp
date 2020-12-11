#include <cmath>
#include <string>
#include <vector>
#include <iostream>

#include "data.hpp"
#include "exceptions.hpp"
#include "matrix.hpp"
#include "vectorf.hpp"
#include "qrfactorization.hpp"

// To fit the tail of the Taylor series, ignore the first kstart terms
// of the TCs.
//
// If Taylor coefficients TC(i) are unscaled, then the scaled TCs are
// T(i) = (1/scale)^i TC(i).
//
// Multiply Rc by scale to get Rc when computing with scaled TCs.
//
// See this from ratio test as Rc is 1/L where
//
//   L = (1/scale)*limit of abs(TC(i+1)/TC(i)) = limit of abs(T(i+1)/T(i)).
//
// To justify the formula for Rc, consider the problem of fitting m and b
//
// log10(abs(T(i))) = y(i) = b + m*i = [1 i](b,m)
//
// by linear least squares. Then
//
// Log10(L) = Log10(abs(T(i+1))/abs(T(i))) =
// Log10(abs(T(i+1)))-Log10(abs(T(i)))
//          = m*(i+1)+b - (m*i+b) = m
//
// It follows that
//
// rc = 1/L = scale/pow(10, m).
int topline(const std::vector<double> &coeff, const double &scale, double &rc,
                 double &order)
{
  int m{(int)coeff.size()-TOPLINE_KSTART}, n{2};
  std::string message;
  matrix<double> W(m,n);
  vectorf<double> beta(n);
  vectorf<double> b(m);

  // Less than 10 coefficients to estimate Rc
  if (m < 10)
  {
    message =
        "Need at least 10 Taylor coefficients to estimate Radius of Convergence\n";
    throw sayMessage(message);
  }

  // Construct the least squares linear system
  for (int k{TOPLINE_KSTART}; k < (int)coeff.size(); k++)
  {
    int row = k - TOPLINE_KSTART + 1;
    W(row, 1) = 1.0;
    W(row, 2) = (double)k;
    b(row) = log10(fabs(coeff[k]));
  }

  // Solve W beta = b for beta
  qr(m, n, W, b, beta);

  rc = scale / pow(10, beta(2));
  order = 1.0;
  // Return 2-norm of residuals
  return 0;
}
