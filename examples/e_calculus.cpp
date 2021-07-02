#include <cmath>
#include <vector>

#include <iostream>
#include "matrix.hpp"
#include "vectorf.hpp"

#include "roc.hpp"

#define NUM_COEFF 30

// Problem from
// https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx
//
// f(x) = sum_0^{infty} \frac{2^n}{n} (4 x -8)^n on 15/8 <= x < 17/8 has Rc=1/8.
//
// We want to fit the tail of the Taylor series, ignore the first 17 terms of
// the TCs.
int main()
{
  // Variables to interact with roc
  std::vector<double> coeffs(NUM_COEFF);
  double scale{1.0};
  double rc{0};
  double order{0};

  for (int k = 0; k < NUM_COEFF; k++)
  {
    coeffs[ k ] = std::pow(8, k + 1) / ((double)(k + 1));
  }
  vectorf<double> ts(coeffs);
  std::cout << "COEFFS =\n" << ts << '\n';
  std::cout << "Scale = " << scale << "\n";
  std::cout << "Size  = " << coeffs.size() << "\n";
  // Call ROC here
  roc(coeffs, scale, rc, order);
}
