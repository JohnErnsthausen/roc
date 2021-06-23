#include <cmath>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "matrix.hpp"
#include "vectorf.hpp"
#include "linearalgebra.hpp"

void constructLinearLeastSquaresRow(const std::vector<double> &coeff, const int k,
                                    double &w1, double &w2, double &b)
{
  w1 = 1.0;
  w2 = (double)k;
  b = std::log10(fabs(coeff[k]));
}

void constructLinearLeastSquaresSystem(const std::vector<double> &coeff,
                                       const int kstart, matrix<double> &W,
                                       vectorf<double> &b)
{
  // Storage required to construct Top-Line system
  if ((int)W.get_rows() != (int)coeff.size() - kstart)
  {
    std::string message{
        "Insufficient storage to construct Top-Line system. Have [" +
        std::to_string(W.get_rows()) + "] Need [" +
        std::to_string((int)coeff.size() - kstart) + "]\n"};
    throw sayMessage(message);
  }
  // Compatible storage
  if( b.get_size() != W.get_rows() )      
  {
    std::string message{
        "Incompatible storage between Matrix and its range in construction of Linear Least Squares system. Have [" +
        std::to_string(W.get_rows()) + "] matrix rows.\nHave [" +
        std::to_string(b.get_size()) + "] rhs rows\n"
    };
    throw sayMessage(message);
  }

  for (int k{kstart}, row{1}; k < (int)coeff.size(); k++, row++)
  {
    constructLinearLeastSquaresRow(coeff, k, W(row, 1), W(row, 2), b(row));
  }
}

// Relative error in all equations starting from TOPLINE_KSTART
double errorTopLine(const std::vector<double> &coeff, const vectorf<double> &beta)
{
  int m = (int)coeff.size() - TOPLINE_KSTART;
  double w1, w2, b;
  vectorf<double> rhs(m);
  vectorf<double> residual(m);
  for (int k{TOPLINE_KSTART}, row{1}; k < (int)coeff.size(); k++, row++)
  {
    constructLinearLeastSquaresRow(coeff, k, w1, w2, b);
    residual(row) = w1 * beta(1) + w2 * beta(2) - b;
    rhs(row) = b;
  }
  double bnrm2 = norm2(m, rhs.data());
  double rnrm2 = norm2(m, residual.data());
  double error = rnrm2/bnrm2;
  return error;
}

// To fit the tail of the Taylor series, ignore the first kstart terms
// of the TCs.
//
// Suppose Taylor coefficients TC(i) are unscaled. The scaled TCs are
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
  int nUse = TOPLINE_NUSE;
  int kStart = TOPLINE_KSTART;
  int m{(int)coeff.size() - kStart}, n{2};
  matrix<double> W(m, n);
  vectorf<double> beta(m);

  // coeff has less than kStart+nUse coefficients to estimate Rc
  if ((int)coeff.size() < kStart+nUse)
  {
    std::string message{
        "Need at least [" + std::to_string(kStart+nUse) +
        "] Taylor coefficients to estimate Radius of Convergence with Top Line Analysis\n"};
    throw sayMessage(message);
  }

  // Construct the least squares linear system
  constructLinearLeastSquaresSystem(coeff, kStart, W, beta);

  // Solve W beta = b for beta
  MinNormSolution(m, n, W.data(), beta.data());

  rc = scale / std::pow(10, beta(2));
  order = std::numeric_limits<double>::quiet_NaN();

  // Compute relative error
  double error = errorTopLine(coeff, beta);
  return error;
}
