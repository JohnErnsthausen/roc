#include "dist.h"
#include <tgmath.h>
#include "mathext.h"
#include "swap.h"

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
  const double zer = 0.0, one = 1.0;

  int ix = 1, iy = 1;
  double diff = 0.0, dx = 0.0, dy = 0.0;
  double scale = 0.0, sum = 0.0, tmp = 0.0, trm = 0.0;

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
    ix = ix + incx;
    if (kvec == 1)  // One vector case
    {
      diff = dx;
    }
    else  // Two vector case
    {
      dy = y(iy);
      iy = iy + incy;

      // Compute diff = dx - dy
      diff = diff_avoids_subtractive_cancellation(dx, dy);
    }

    // Sum the scaled squares (all less than or equal to one) of the nonzero
    // terms.

    if (diff != zer)
    {
      trm = fabs(diff);
      if (scale < trm)
      {
        tmp = scale / trm;
        sum = one + sum * tmp * tmp;
        scale = trm;
      }
      else
      {
        tmp = trm / scale;
        sum += tmp * tmp;
      }
    }
  }
  return scale * sqrt(sum);
}

// If dx and dy are not of the same sgn, then safely subtract
//
// else
//
// If dx and dy are of the same sign, then subtract by
// 1.) make sure dx > dy
// 2.a) We have sgn(dx)=sgn(dy).
//      If dx=0, then dy=0. Thus diff=0.
// 2.b) else Subtract dx - dy = dx*(one-dy/dx)
double diff_avoids_subtractive_cancellation(double dx, double dy)
{
  const double zer = 0.0, one = 1.0;
  double diff = 0.0;

  if (sgn(dx) != sgn(dy))
  {
    diff = dx - dy;
  }
  else
  {
    if (fabs(dx) < fabs(dy)) dswap(&dx, &dy);
    if (dx == zer)
    {
      diff = zer;
    }
    else
    {
      diff = dx * (one - dy / dx);
    }
  }
  return diff;
}
