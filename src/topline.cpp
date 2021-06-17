#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "matrix.hpp"
#include "vectorf.hpp"
#include "linearalgebra.hpp"

void constructLinearLeastSquaresSystem(const std::vector<double> &coeff,
                                       const int kstart, matrix<double> &W,
                                       vectorf<double> &b)
{
  // Storage required to construct Top-Line system
  if ((int)coeff.size() - kstart > (int)W.get_rows())
  {
    std::string message{
        "Insufficient storage to construct Top-Line system. Have [" +
        std::to_string(W.get_rows()) + "] Need [" +
        std::to_string((int)coeff.size() - kstart) + "]\n"};
    throw sayMessage(message);
  }

  for (int k{kstart}; k < (int)coeff.size(); k++)
  {
    int row = k - kstart + 1;
    W(row, 1) = 1.0;
    W(row, 2) = (double)k;
    b(row) = log10(fabs(coeff[ k ]));
  }
}

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
double topline(const std::vector<double> &coeff, const double &scale,
               double &rc, double &order)
{
  int m{(int)coeff.size() - TOPLINE_KSTART}, n{2};
  matrix<double> W(m, n);
  vectorf<double> beta(n);
  vectorf<double> b(m);

  // Less than TOPLINE_NUSE coefficients to estimate Rc
  if (m < TOPLINE_NUSE)
  {
    std::string message{
        "Need at least [" + std::to_string(TOPLINE_NUSE) +
        "] Taylor coefficients to estimate Radius of Convergence\n"};
    throw sayMessage(message);
  }

  // Construct the least squares linear system
  constructLinearLeastSquaresSystem(coeff, TOPLINE_KSTART, W, b);

  // Save a copy of linear system for computing residual
  matrix<double> WSaved{W};
  vectorf<double> bSaved{b};

  // Solve W beta = b for beta
  MinNormSolution(m, n, W.data(), beta.data());

  rc = scale / pow(10, beta(2));
  order = std::numeric_limits<double>::quiet_NaN();

  // Compute and return 2-norm of residuals
  for (int i{1}; i <= m; i++)
    bSaved(i) -= WSaved(i, 1) * beta(1) + WSaved(i, 2) * beta(2);
  return norm2(m, bSaved.data());
}
