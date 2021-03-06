//  gcc -I.. GetTS.cpp -o ex.exe -lstdc++ -lm
//  ./ex.exe
//
#include <cmath>
#include <iostream>
#include "tadiff.h"
using namespace std;
using namespace fadbad;

T<double> func(const T<double>& x)
{
  T<double> y = 1.0 + 25.0 * x * x;
  y = pow(y, 3.1415);

  // T<double> y = 1.0 - x;
  return 1.0 / y;
}

int main()
{
  double x0{0.3}, hbar{1.0};

  T<double> x, f;
  x = x0;         // Expansion point
  x[ 1 ] = hbar;  // Taylor-expand wrt. x (dx/dx=1)
  f = func(x);    // Evaluate function and record DAG
  f.eval(30);     // f[0]...f[30] now contains the Taylor-coefficients.
  cout.precision(16);
  cout << scientific;
  for (int i = 0; i < 30; i++)
  {
    cout << "    " << f[ i ] << '\n';  // The i'th taylor coefficient
  }
  return 0;
}
