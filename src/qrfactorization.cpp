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
  vector<int> ipiv(n, 0);
  vector<double> tau(n, 0.0);
  vector<double> wrk(m, 0.0);
  double safmin{DBL_MIN};

  // Factor
  qrf(m, n, &W(1, 1), m, &ipiv[ 0 ], &tau[ 0 ], &wrk[ 0 ], safmin, &ier);
  if (ier != 0)
  {
    message = "QRFactorization error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
    return ier;
  }
  // Solve
  qrs(m, n, &W(1, 1), m, &tau[ 0 ], &b(1), &x(1), &ier);
  if (ier != 0)
  {
    message = "QRSolver error with ier= " + to_string(ier) + " \n";
    throw sayMessage(message);
    return ier;
  }
  // Multiply solution x by permutation and store it in x
  for (int i{0}; i < n; i++)
  {
    int j = ipiv[ i ];
    tau[ i ] = x(j);
  }
  for (int i{0}; i < n; i++)
  {
    int j = i + 1;
    x(j) = tau[ i ];
  }
  return ier;
}
