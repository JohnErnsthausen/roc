#include <cmath>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "matrix.hpp"
#include "qrfactorization.hpp"
#include "vectorf.hpp"

void constructSixTermSystem(const vectorf<double> &coeff, const int nUse,
                            matrix<double> &W, vectorf<double> &b)
{
  // TODO Check coeff.get_size, nUse, and W.get_rows are compatible?

  for (int i{1}, k{(int)coeff.get_size() - nUse}; i <= nUse; i++, k++)
  {
    W(i, 1) = 2.0 * coeff(k);
    W(i, 2) = 2.0 * (k - 1) * coeff(k);
    W(i, 3) = -2.0 * coeff(k - 1);
    W(i, 4) = -(k - 2) * coeff(k - 1);
    b(i) = k * coeff(k + 1);
  }
}

void testBeta4(double beta4)
{
  if (beta4 < 0)
  {
    std::string message =
        "Unconstrained optimization lead to Sqrt of negative number: " +
        std::to_string(beta4) + " \n";
    throw sayMessage(message);
  }
}

void testRCSix(double rc)
{
  if (std::isnan(rc))
  {
    std::string message =
        "The radius of convergence is infinity, which is highly unlikely\n";
    throw sayMessage(message);
  }
}

void testCosTheta(double cosTheta)
{
  if (std::isnan(cosTheta))
  {
    std::string message =
        "Unconstrained optimization lead to infinite CosTheta which is not in "
        "[-1, 1]\n";
    throw sayMessage(message);
  }
  if ((cosTheta < -1.0) || (cosTheta > 1.0))
  {
    std::string message = "Unconstrained optimization lead to CosTheta [" +
                          std::to_string(cosTheta) + "] not in [-1, 1]\n";
    throw sayMessage(message);
  }
}

double testSingularityOrder(double singularityOrder1, double singularityOrder2)
{
  double order{0.0};
  if (std::isnan(singularityOrder1) && std::isnan(singularityOrder2))
  {
    std::string message =
        "Unconstrained optimization lead to NaN for Order of Singularity\n";
    throw sayMessage(message);
  }

  if (std::isnan(singularityOrder1) && !std::isnan(singularityOrder2))
    order = singularityOrder2;

  if (!std::isnan(singularityOrder1) && std::isnan(singularityOrder2))
    order = singularityOrder1;

  if (!std::isnan(singularityOrder1) && !std::isnan(singularityOrder2))
    order = (singularityOrder1 + singularityOrder2) / 2.0;

  return order;
}

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
double sixterm(const std::vector<double> &coeff, const double &scale,
               double &rc, double &order)
{
  int nUse = SIXTERM_NUSE;
  int m{nUse}, n{4};
  matrix<double> W(m, n);
  vectorf<double> beta(n);
  vectorf<double> b(m);
  vectorf<double> tc(coeff);

  constructSixTermSystem(tc, nUse, W, b);

  // Solve W beta = b for beta
  qr(m, n, W, b, beta);

  // Interpret the variables found from Least Squares Optimization Problem
  double hOverRc{0.0}, cosTheta{0.0}, singularityOrder1{0.0},
      singularityOrder2{0.0};

  // Evaluate h/Rc
  testBeta4(beta(4));
  hOverRc = sqrt(beta(4));

  // Evaluate Rc
  rc = scale / hOverRc;
  testRCSix(rc);

  // Evaluate Cos(Theta)
  cosTheta = beta(2) / hOverRc;
  testCosTheta(cosTheta);

  // Evaluate Order of the Singularity
  singularityOrder1 = beta(1) / beta(2);
  singularityOrder2 = beta(3) / beta(4);
  order = testSingularityOrder(singularityOrder1, singularityOrder2);

  // Compare order
  // printf( "Abs of difference between two computations for order:
  // [%22.16f]\n",
  //     fabs(singularityOrder1-singularityOrder2));
  // TODO Add this error to checked error on return?

  // Evaluate previous W equation at the solution of the least squares
  // optimization problem. Return the backward error, the absolute value of this
  // evaluation.
  //
  // TODO Use all equations?
  nUse++;
  int k = coeff.size() - nUse;
  double check =
      k * tc(k + 1) -
      ((2.0 * tc(k)) * beta(1) + (2.0 * (k - 1) * tc(k)) * beta(2) +
       (-2.0 * tc(k - 1)) * beta(3) + (-(k - 2) * tc(k - 1)) * beta(4));

  return fabs(check);
}
