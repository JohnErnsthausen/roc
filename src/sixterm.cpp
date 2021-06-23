#include <cmath>
#include <string>
#include <vector>

#include <iostream>

#include "data.hpp"
#include "exceptions.hpp"
#include "matrix.hpp"
#include "vectorf.hpp"
#include "linearalgebra.hpp"

void constructSixTermRow(const vectorf<double> &coeff, const int k,
                         double &w1, double &w2, double &w3, double &w4, double &b)
{
  w1 = 2.0 * coeff(k);
  w2 = 2.0 * (k - 1) * coeff(k);
  w3 =-2.0 * coeff(k - 1);
  w4 =-(k - 2) * coeff(k - 1);
  b  = k * coeff(k + 1);
}

void constructSixTermSystem(const vectorf<double> &coeff, const int from, const int to,
                            matrix<double> &W, vectorf<double> &b)
{
  // Suggicient storage. Remember it must be a contiguous block.
  if( (to-from+1 != (int)W.get_rows()) || ((int)W.get_cols() != 4) )      
  {
    std::string message{
        "Insufficient storage to construct Six Term system. Have [" +
        std::to_string(W.get_rows()) + "] rows Need [" +
        std::to_string(to-from+1) + "] rows\nHave [" +
        std::to_string(W.get_cols()) + "] cols Need [4] cols\n"
    };
    throw sayMessage(message);
  }
  // Compatible storage
  if( b.get_size() != W.get_rows() )      
  {
    std::string message{
        "Incompatible storage between Matrix and its range in construction of Six Term system. Have [" +
        std::to_string(W.get_rows()) + "] matrix rows.\nHave [" +
        std::to_string(b.get_size()) + "] rhs rows\n"
    };
    throw sayMessage(message);
  }
  for (int i{1}, k{from}; k <= to; i++, k++)
  {
    constructSixTermRow(coeff, k, W(i, 1), W(i, 2), W(i, 3), W(i, 4), b(i));
  }
  // std::cout << "Rows [" + std::to_string(W.get_rows()) + "]\n";
  // std::cout << "Cols [" + std::to_string(W.get_cols()) + "]\n";
  // std::cout << "W =\n" << W << '\n';
  // std::cout << "b =\n" << b << '\n';
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

  double diff = fabs(singularityOrder1 - singularityOrder2);
  if ( diff > TOL )
  {
    std::cout.precision(16);
    std::cout << std::scientific << "|s1 - s2| = " + std::to_string(diff) + "\n";
  }

  return order;
}

// Relative error in all equations starting from SIXTERM_KSTART
double errorSixTerm(const vectorf<double> &coeff, const vectorf<double> &beta)
{
  int dim{0};
  double w1, w2, w3, w4;
  vectorf<double> residual(coeff.get_size());
  vectorf<double> rhs(coeff.get_size());
  for (int i{1}, k{SIXTERM_KSTART}; k < (int)coeff.get_size(); i++, k++)
  {
    dim = i;
    constructSixTermRow(coeff, k, w1, w2, w3, w4, rhs(i));
    residual(i) = w1*beta(1) + w2*beta(2) + w3*beta(3) + w4*beta(4) - rhs(i);
  }
  double bnrm2 = norm2(dim, rhs.data());
  double rnrm2 = norm2(dim, residual.data());
  double error = rnrm2/bnrm2;
  return error;
}

// The six-term-test of Chang and Corliss
//
// The Six-Term model subproblem is the Six-Term model at indicies FROM to TO
// where FROM = coeff.size()-nUse and TO = coeff.size()-1 and nUse is defined in
// header file data.h
//
// The method returns the relative error of the Six-Term model evaluated
// at the computed solution of the Six-Term model subproblem from index
// SIXTERM_KSTART to index coeff.size()-1.
double sixterm(const std::vector<double> &coeff, const double &scale,
               double &rc, double &order)
{
  int nUse = SIXTERM_NUSE;
  int kStart = SIXTERM_KSTART;
  int m{nUse}, n{4};
  matrix<double> W(m, n);
  vectorf<double> beta(m);
  vectorf<double> tc(coeff);

  // coeff has less than kStart+nUse coefficients to estimate Rc
  if ((int)coeff.size() < kStart+nUse)
  {
    std::string message{
        "Need at least [" + std::to_string(kStart+nUse) +
        "] Taylor coefficients to estimate Radius of Convergence with Six Term Analysis\n"};
    throw sayMessage(message);
  }

  // from must be greater than 0
  int from = (int)coeff.size()-nUse;
  // to must be less than coeff.size()
  int to = (int)coeff.size()-1;
  constructSixTermSystem(tc, from, to, W, beta);

  // Solve W beta = b for beta
  MinNormSolution(m, n, W.data(), beta.data());

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

  return errorSixTerm(tc, beta);
}

