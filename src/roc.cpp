#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "exceptions.hpp"
#include "roc.hpp"

using namespace std;

// The vector coeff is expected to have at least length MINTERMS (put to
// data.hpp)
int roc(const vector<double> &coeff, const double &scale, double &rc,
        double &order)
{
  string message;
  double hrc{0.0}, hrc_check{0.0};

  rc = 0.0;
  rc *= scale;
  order = 0.0;

  // Check for sufficient data to perform three term analysis
  int k = (int)coeff.size();
  if (k < MINTERMS)
  {
    message = "The number of coefficients [" + to_string((int)coeff.size()) +
              "] must be greater than [" + to_string(MINTERMS) + "].\n";
    throw sayMessage(message);
  }

  // Check for agreement between computations against TOL
  //
  // If no agreement, three term analysis is said to have failed
  // because the coefficients do not represent a pole.
  if (fabs(hrc - hrc_check) > TOL)
  {
    cout.precision(16);
    cout << scientific;
    cout << "|hrc - hrc_check| = " << fabs(hrc - hrc_check) << "\n";
  }

  return 1;
}
