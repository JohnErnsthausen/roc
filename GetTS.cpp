//  gcc -I.. GetTS.cpp -o ex.exe -lstdc++ -lm
//  ./ex.exe
//
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>

#include "tadiff.h"

#define DIM 35

#define prec double

using namespace fadbad;

T<prec> func(const T<prec>& x)
{
  T<prec> z = 1. + 25.0 * x * x;
  z = 1.0/z;
  //z = 1.0/(z*z*z*z*z*z*z*z);
  //return z;

  // x[0] = 1.0/(1.0+25.0*t*t) + 1.0/(1.0+25.0*t*t)^2;
  std::vector<T<prec>> powers(DIM, 1.0);
  for(int i{0}; i<DIM; i++)
  {
    for(int j{0}; j<=i; j++)
    {
      powers[j] *= z;
    }
  }
  reverse(powers.begin(), powers.end());
  return std::accumulate( powers.begin(), powers.begin()+8, T<prec> (0.0) );
}

int main()
{
  double x0{-0.3}, hbar{1.0};

  T<prec> x, f;
  x = x0;         // Expansion point
  x[ 1 ] = hbar;  // Taylor-expand wrt. x (dx/dx=1)
  f = func(x);    // Evaluate function and record DAG
  f.eval(DIM);     // f[0]...f[30] now contains the Taylor-coefficients.
  std::cout.precision(16);
  std::cout << std::scientific;
  for (int i{0}; i < DIM; i++)
  {
    std::cout << "    " << (double)f[ i ] << '\n';  // The i'th taylor coefficient
  }
  return 0;
}
