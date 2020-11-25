#include <cmath>
#include <iostream>
#include <vector>

#include "exceptions.hpp"
#include "roc.hpp"

using namespace std;

// The vector coeff is expected to have at least length MINTERMS (put to
// data.hpp)
int roc(const vector<double> &coeff, const double &scale, double &rc,
        double &order)
{
  double hrc{0.0}, hrc_check{0.0};

  rc = 0.0;
  rc *= scale;
  order = 0.0;

  // Check for sufficient data to perform three term analysis
  int k = (int)coeff.size();
  if (k < MINTERMS) throw morecoefficients();

  // Check for agreement between computations against TOL
  //
  // If no agreement, three term analysis is said to have failed
  // becuase the coefficients do not represent a pole.
  if (fabs(hrc - hrc_check) > TOL)
  {
    cout.precision(16);
    cout << scientific;
    cout << "|hrc - hrc_check| = " << fabs(hrc - hrc_check) << "\n";
  }

  return 1;
}
