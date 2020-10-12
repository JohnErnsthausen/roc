#include <tgmath.h>
#include <stdio.h>

#include "roc.h"

#define NUM_TC 30

// Problem from
// https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx
//
// f(x) = sum_0^{infty} \frac{2^n}{n} (4 x -8)^n on 15/8 <= x < 17/8 has Rc=1/8.
//
// We want to fit the tail of the Taylor series, ignore the first 17 terms of
// the TCs.
int main(void)
{
  int num_tc = NUM_TC, kstart = 18;
  int res;
  double scale = 1.0;
  double rc = 0.0, slope=0.0, intercept=0.0;
  double tc[NUM_TC] = {0};

  for(int k=0; k<num_tc; k++) { tc[k] = pow(8,k+1)/((double) (k+1)); }

  res = roc(num_tc, tc, scale, kstart, &rc, &slope, &intercept);
  if( res!= 0 )
  {
    printf("Error finding the Radius of Convergence\n");
  }

  printf("The Radius of Convergence Rc for the given data and scale\nRc[%22.16f] scale[%22.16f]\n",
      rc, scale);
  printf("The Least Squares best fit is y = %22.16f x +  %22.16f\n",slope,intercept);
}
