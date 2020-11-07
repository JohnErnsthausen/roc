#include "dist.h"
#include <tgmath.h>
#include "mathext.h"
#define swap(T,x,y) {T tmp = (x); x = (y); y = tmp;}

#define map1(i) (i) - 1
#define x(i) x[ map1(i) ]
#define y(i) y[ map1(i) ]

// Computes either the Euclidean distance between two N-dimensional
// vectors X and Y or the Euclidean norm of one such vector X.
//
// Call the routine
// either
//    with KVEC = 2 and two vectors X and Y of dimension N stored
//    with storage-increments INCX and INCY, respectively.
// or
//    with KVEC = 1 and one vector X of dimension N stored
//    with storage-increment INCX. In this case the second array
//    is not referenced and can be a dummy array or simply the
//    array X again.
//
// If N .LE. 0 then zero is returned, if N .GE. 1 then the
// storage increments cannot be zero.
//
// This code is inspired by the LAPACK function DNRM2 written by
// Sven Hammarling.
//
// Variables in the calling sequence
// ---------------------------------
//    N     I    IN   Dimension of the vectors X and Y
//    X     D    IN   The first vector of dimension N
//    INCX  I    IN   Storage increment of X
//    Y     D    IN   The second vector of dimension N
//    INYY  D    IN   Storage increment of Y
//    KVEC  I    IN   Number of vectors
//                    KVEC = ANY INTEGER OTHER THAN 2
//                              Only one vector, namely X, is given,
//                              the Euclidean norm of X is computed
//                              and the Y array is not referenced
//                              This is the default
//                    KVEC = 2  Two vectors X and Y  are given,
//                              the Euclidean distance between
//                              X and Y is computed
//
double ddist2(int n, double *x, int incx, double *y, int incy, int kvec)
{
  const double zer = 0.0;
  int ix = 1, iy = 1;
  double diff = 0.0, dx = 0.0, dy = 0.0;
  double scale = 0.0, sum = 0.0;

  if (kvec != 2) kvec = 1;

  // Define the inner product to be zero whenever the dimension is not positive
  // or a storage incrementor (incx/incy) is zero
  if (n <= 0) return zer;
  if (kvec == 1 && incx == 0) return zer;
  if (kvec == 2 && incx == 0 && incy == 0) return zer;

  // Initialize incrementors

  ix = 1;
  if (incx < 0) ix = (-n + 1) * incx + 1;
  iy = 1;
  if (incy < 0) iy = (-n + 1) * incy + 1;

  sum = zer;
  scale = zer;
  for (int j = 1; j <= n; j++)
  {
    dx = x(ix);
    ix += incx;
    if (kvec == 1)  // One vector case
    {
      diff = dx;
    }
    else  // Two vector case
    {
      dy = y(iy);
      iy += incy;
      diff = dx - dy;
    }

    // Scaled squares (less than or equal to one) of the nonzero terms.
    add_next_element_squared(diff, &sum, &scale);
  }
  return scale * sqrt(sum);
}

// Computes the Euclidean norm of X.
//
// If N .LE. 0 then zero is returned, if N .GE. 1 then the
// storage increments cannot be zero.
//
// This code is inspired by the LAPACK function DNRM2 written by
// Sven Hammarling.
//
// Variables in the calling sequence
// ---------------------------------
//    N     I    IN   Dimension of the vector X
//    X     D    IN   Vector of dimension N
//    INCX  I    IN   Storage spacing between elements of X
//
double dnrm2(int n, double *x, int incx)
{
  const double zer = 0.0;
  int ix = 1;
  double dx = 0.0;
  double scale = 0.0, sum = 0.0;

  // Define the norm2 to be zero whenever the dimension is not positive
  // or a storage incrementor is zero
  if (n < 1 || incx==0) return zer;

  // Initialize incrementors
  if (incx < 0) ix = (-n + 1) * incx + 1;

  for (int j = 1; j <= n; j++)
  {
    dx = x(ix);
    ix += incx;
    // Scaled squares (less than or equal to one) of the nonzero terms.
    add_next_element_squared(dx, &sum, &scale);
  }
  return scale * sqrt(sum);
}

void add_next_element_squared(double xi, double *sum, double *scale)
{
  double bar = *scale;
  double ssq = *sum;
  const double zer = 0.0, one = 1.0;
  double tmp = 0.0, absxi = 0.0;

  // If xi==0, then the remaining statements do not permit a division by zero
  // because bar >= 0 implies that bar < 0 is always false.
  // However a quick return prohibits a wasted multiplication and sum.
  if (xi == zer) return;
  
  absxi = fabs(xi);
  if (bar < absxi)
  {
    tmp = bar / absxi;
    ssq = one + ssq * tmp * tmp;
    bar = absxi;
  }
  else
  {
    tmp = absxi / bar;
    ssq += tmp * tmp;
  }

  *scale = bar;
  *sum = ssq;
}
