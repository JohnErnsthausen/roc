#include <cfloat>
#include <iostream>
#include <string>
#include <vector>

extern "C"
{
#include "qrfactorization.h"
}
#include "exceptions.hpp"
#include "matrix.hpp"
#include "qrfactorization.hpp"
#include "vectorf.hpp"

using namespace std;

int qr(const int m, const int n, matrix<double> &W, vectorf<double> &b,
       vectorf<double> &x)
{
  int ier{0};
  string message;
  vectorf<int> ipiv(n);
  vectorf<double> tau(n);

  // Factor
  factor(m, n, W, tau, ipiv);
  // Solve
  solve(m, n, W, tau, b, x);
  // Multiply solution x by permutation and store it in x
  permute(n, x, ipiv, tau);

  return ier;
}

int factor(const int m, const int n, matrix<double> &W, vectorf<double> &tau,
           vectorf<int> &ipiv)
{
  int ier{0};
  double safmin{DBL_MIN};
  vectorf<double> wrk(m);
  qrf(m, n, &W(1, 1), m, &ipiv(1), &tau(1), &wrk(1), safmin, &ier);
  if (ier != 0)
  {
    std::string message =
        "This QRFactorization error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
  }
  return ier;
}

int solve(const int m, const int n, matrix<double> &W, vectorf<double> &tau,
          vectorf<double> &b, vectorf<double> &x)
{
  int ier{0};
  qrs(m, n, &W(1, 1), m, &tau(1), &b(1), &x(1), &ier);
  if (ier != 0)
  {
    std::string message = "QRSolver error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
  }
  return ier;
}

int permute(const int n, vectorf<double> &x, vectorf<int> &ipiv,
            vectorf<double> &wrk)
{
  int ier{0};
  for (int i{1}; i <= n; i++)
  {
    wrk(i) = x(ipiv(i));
  }
  for (int i{1}; i <= n; i++)
  {
    x(i) = wrk(i);
  }
  return ier;
}
