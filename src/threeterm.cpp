#include <cmath>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "matrix.hpp"
#include "vectorf.hpp"
#include "linearalgebra.hpp"

void constructThreeTermRow(const vectorf<double> &coeff, const int k,
                           double &w1, double &w2, double &b)
{
  w1 = (k - 1) * coeff(k);
  w2 = coeff(k);
  b  = k * coeff(k + 1);
}

void constructThreeTermSystem(const vectorf<double> &coeff, const int from, const int to,
                            matrix<double> &W, vectorf<double> &b)
{
  for (int i{1}, k{from}; k <= to; i++, k++)
  {
    constructThreeTermRow(coeff, k, W(i, 1), W(i, 2), b(i));
  }
  // std::cout << "Rows [" + std::to_string(W.get_rows()) + "]\n";
  // std::cout << "Cols [" + std::to_string(W.get_cols()) + "]\n";
  // std::cout << "W =\n" << W << '\n';
  // std::cout << "b =\n" << b << '\n';
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

// Relative error in all equations starting from THREETERM_KSTART
double errorThreeTerm(const vectorf<double> &coeff, const vectorf<double> &beta)
{
  int dim{0};
  double w1, w2;
  vectorf<double> residual(coeff.get_size());
  vectorf<double> rhs(coeff.get_size());
  for (int i{1}, k{THREETERM_KSTART}; k < (int)coeff.get_size(); i++, k++)
  {
    dim = i;
    constructThreeTermRow(coeff, k, w1, w2, rhs(i));
    residual(i) = w1*beta(1) + w2*beta(2) - rhs(i);
  }
  double bnrm2 = norm2(dim, rhs.data());
  double rnrm2 = norm2(dim, residual.data());
  double error = rnrm2/bnrm2;
  return error;
}

// The three-term-test of Chang and Corliss
double threeterm(const std::vector<double> &coeff, const double &scale,
                 double &rc, double &order)
{
  int nUse = THREETERM_NUSE;
  int m{nUse}, n{2};
  matrix<double> W(m, n);
  vectorf<double> beta(m);
  vectorf<double> tc(coeff);

  // from must be greater than 0
  int from = (int)coeff.size()-nUse;
  // to must be less than coeff.size()
  int to = (int)coeff.size()-1;
  constructThreeTermSystem(tc, from, to, W, beta);

  // Solve W beta = b for beta
  MinNormSolution(m, n, W.data(), beta.data());

  // Interpret the variables found from Least Squares Optimization Problem
  double hOverRc = beta(1);
  rc = scale / hOverRc;
  testRCThree(rc);

  order = beta(2) / beta(1);
  testOrder(order);

  return errorThreeTerm(tc, beta);
}
