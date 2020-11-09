#include <stdio.h>
#include <tgmath.h>

#include "dist.h"
#include "mathext.h"
#include "qrfactorization.h"

#define map(i, j) ((j)-1) * lda + ((i)-1)
#define a(i, j) a[ map(i, j) ]
#define map1(i) (i) - 1
#define x(i) x[ map1(i) ]
#define y(i) y[ map1(i) ]
#define ipiv(i) ipiv[ map1(i) ]
#define wrk(i) wrk[ map1(i) ]
#define tau(i) tau[ map1(i) ]
#define swap(T, x, y) \
  {                   \
    T tmp = (x);      \
    x = (y);          \
    y = tmp;          \
  }

// Compute the QR factorization with column pivoting on all columns
// of an M by N matrix A
//                     A*P = Q*R
//
// The matrix Q is represented as a product of elementary reflectors
//
//    Q = H(1) H(2) . . . H(k),   where k = min(m,n).
//
// Each H(i) has the form
//
//    H(i) = I - tau * v * v'
//
// where tau is a real scalar, and v is a real vector with v(1:i-1) = 0
// and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
// and tau in TAU(i).
//
// The matrix P is represented in IPIV as follows: If IPIV(j) = i
// then the jth column of P is the ith canonical unit vector.
//
// Variables in the calling sequence
// ---------------------------------
// M      I   IN   The number of rows of the matrix A.  M >= 0.
// N      I   IN   The number of columns of the matrix A.  N >= 0.
// A      D   IN   The given M by N matrix
//            OUT  The elements on and above the diagonal contain
//                 the min(m,n) by n upper trapezoidal matrix R
//                 (R is upper triangular if m >= n); the elements
//                  belowthe diagonal, together with the array TAU,
//                  represent the orthogonal matrix Q as a product
//                  of min(m,n) elementary reflectors
// LDA    I   IN   The leading dimension of A. LDA >= max(1,M).
// IPIV   I   OUT  If IPIV(i) = k, then the i-th column of A*P was the
//                  k-th column of A.
// TAU    I   OUT  Array of dimension min(M,N), The scalar factors of
//                 the elementary reflectors
// WRK    D   WK   Work array of  dimension N
// SAFMIN D   IN   Safe minimum such that 1.0/SAFMIN does not
//                 overflow
// IER    I   OUT  Error indicator
//                 IER = 0   successful exit
//                 IER = 1   data-error
//                 IER = 2   input data error in HOUSL
//                 IER < 0   zero pivot encountered
int qrf(int m, int n, double *a, int lda, int *ipiv, double *tau, double *wrk,
        double safmin, int *ier)
{
  const double fact = 0.05, one = 1.0, zer = 0.0;
  int imax = 0, ip1 = 0, mn = 0;
  double cnrm = 0.0, cnrmj = 0.0, taui = 0.0, tmp1 = 0.0, tmp2 = 0.0,
         tmp3 = 0.0;

  // Validate input
  if ((m < 0) || (n < 0) || (lda < m))
  {
    *ier = 1;
    return 1;
  }

  // Quick return
  if ((m == 0) || (n == 0))
  {
    *ier = 0;
    return 0;
  }

  // Initialize column norms and pivot array
  for (int i = 1; i <= n; i++)
  {
    cnrm = dnrm2(m, a + map(1, i), 1);
    tau(i) = cnrm;
    wrk(i) = cnrm;
    ipiv(i) = i;
  }

  mn = n;
  if (m < n) mn = m;

  // Main loop for the factorization
  for (int i = 1; i <= mn; i++)
  {
    ip1 = i + 1;
    // determine pivot column
    cnrm = tau(i);
    imax = i;
    if (i < n)
    {
      for (int j = ip1; j <= n; j++)
      {
        if (tau(j) > cnrm)
        {
          imax = j;
          cnrm = tau(j);
        }
      }
      if (cnrm == zer)
      {
        printf("Message from QRF: Zero pivot encountered\n");
        *ier = -i;
      }
      // Swap the columns if necessary
      if (imax != i)
      {
        for (int j = 1; j <= m; j++)
        {
          swap(double, a(j, imax), a(j, i));
        }
        swap(int, ipiv(imax), ipiv(i));
        tau(imax) = tau(i);
        wrk(imax) = wrk(i);
      }
    }
    if (cnrm == zer)
    {
      printf("Message from QRF: Zero pivot encountered\n");
      *ier = -i;
    }
    else
    {
      // Generate elementary reflector H(i)
      taui = zer;
      if (i < m)
        housg(m - i + 1, a + map(i, i), a + map(ip1, i), 1, &taui, safmin);
      if (i < n)
      {
        // Apply h(i) to a(i:m,i+1:n) from the left
        housl(m - i + 1, n - i, a + map(i, i), 1, taui, a + map(i, ip1), lda,
              ier);
        if (*ier != 0)
        {
          printf("Message from QRF: Input-dimension error in housl\n");
          *ier = 2;
          return *ier;
        }
        // Update column norms
        for (int j = ip1; j <= n; j++)
        {
          cnrmj = tau(j);
          if (cnrmj != zer)
          {
            tmp1 = fabs(a(i, j)) / cnrmj;
            tmp2 = one - tmp1 * tmp1;
            if (tmp2 <= zer)
            {
              tmp2 = zer;
              tmp3 = zer;
            }
            else
            {
              tmp3 = sqrt(tmp2);
            }
            tmp1 = cnrmj / wrk(j);
            tmp2 = one + fact * tmp2 * tmp1 * tmp1;
            if (tmp2 == one)
            {
              tau(j) = dnrm2(m - i, a + map(ip1, j), 1);
              wrk(j) = tau(j);
            }
            else
            {
              tau(j) = tau(j) * tmp3;
            }
          }
        }
      }
      tau(i) = taui;
    }
  }
  *ier = 0;
  return *ier;
}

// For an  M x N matrix A with M >= N and rank A = N, and a given
// M-vector Y, compute the least squares solution
//
//    min { || A*X - Y || ; X in R^N }
//
// by the algorithm
//
//     1.  Y := Q^T Y
//     2.  IF( M >= N ) THEN Solve R*Z := (I, 0)Y, Set X := Z
//                      ELSE Solve R*Z := Y, Set X = (Z,0)
//
// under the assumption that the arrays A and TAU contain the
// QR-factorization of an M x N matrix A.
//
// The routine returns Q^T Y in the array Y and the least squares
// solution in the array X. The array X may be identified with Y
// if Q^T*Y is not needed.
//
// Variables in the calling sequence
// ---------------------------------
// M    I   IN   Number of rows of the matrix A, M >= N
// N    I   IN   Number of columns of the matrix A
// A    D   IN   Array of dimension LDA x n, the factored matrix
// LDA  I   IN   Leading dimension of A, LDA >= M
// TAU  D   IN   Array of dimension N containing the scalar
//               factor of the elementary reflectors, as
//               returned by QRF
// Y    D   IN   Array of dimension M, the given vector Y
//          OUT  The vector Q^T Y
// X    D   OUT  Array of dimension N, the computed solution
// IER  I   OUT  Error indicator
//               IER = 0   no error
//               IER = 1  input data error
//               IER < 0  zero pivot encountered
//                        IER = -J signifies that the
//                        J-th diagonal element of the
//                        triangular matrix R is zero.
//
int qrs(int m, int n, double *a, int lda, double *tau, double *y, double *x,
        int *ier)
{
  const double zer = 0.0;
  int jm1, jp1;
  double sum, t;

  // Validate input: lda >= m >= 0 holds from these conditions
  if ((m < 0) || (n < 0) || (m < n) || (lda < m))
  {
    *ier = 1;
    return 1;
  }

  // Quick return: 0 = m >= n >= 0 implies m=n=0 so that lda >= 0.
  // If lda=0, then there is a quick return and matrix a is not accessed.
  if ((n == 0) || ((m == 0) && (n == 0)))
  {
    *ier = 0;
    return 0;
  }

  // The case m = 1, under the assumption m >= n, implies n=1
  if (m == 1)
  {
    if (a(1, 1) == zer)
    {
      printf("Message from QRS: Zero pivot encountered\n");
      *ier = -1;
      return -1;
    }
    else
    {
      x(1) = y(1) / a(1, 1);
      *ier = 0;
      return 0;
    }
  }

  // Compute transpose(Q)*y
  //
  // The matrix Q is represented as a product of elementary reflectors
  //
  //    Q = H(1) H(2) . . . H(k),   where k = min(m,n).
  //
  // Each H(i) has the form
  //
  //    H(i) = I - tau * v * v'
  //
  //  v(1) = 1. transpose(Q) = H(j) H(j-1) . . . H(1).
  for (int j = 1; j <= n; j++)
  {
    jp1 = j + 1;

    // Multiply transpose(v)*y
    sum = y(j);  // The first element of the elementary vector is 1
    if (j < m)
    {
      for (int i = jp1; i <= m; i++)
      {
        sum = sum + a(i, j) * y(i);  // v stored in lower triangle of A
      }
    }

    // y overwritten by y - tau(j)*(transpose(v)*y)*v
    if (sum != zer)  // v is not zero
    {
      t = -tau(j) * sum;
      y(j) = y(j) + t;  // The first element of v is 1
      if (j < m)
      {
        for (int i = jp1; i <= m; i++)
        {
          y(i) = y(i) + t * a(i, j);
        }
      }
    }
  }

  // Set x := the first n components of y
  for (int j = 1; j <= n; j++)
  {
    x(j) = y(j);
  }
  // R is stored in the upper triangle of A including diagonal.
  //
  // Solve R*Z := (I,0)Y, and set X := Z
  // or solve  R*Z := Y, and set X := (Z,0)
  for (int j = n; j >= 1; j--)
  {
    if (a(j, j) == zer)
    {
      printf("Message from QRS: Zero pivot encountered\n");
      *ier = -j;
      return -j;
    }
    // Back substitution
    x(j) = x(j) / a(j, j);
    if (j > 1)
    {
      jm1 = j - 1;
      t = -x(j);
      for (int i = 1; i <= jm1; i++)
      {
        x(i) = x(i) + t * a(i, j);
      }
    }
  }
  *ier = 0;
  return *ier;
}

// Generates an n-dimensional Householder reflector
//
//    H = I - tau*( 1 ) * ( 1 v' ),    H' * H = I,   tau scalar
//                ( v )
//
// such that
//
//    H * ( alpha ) = ( beta ),    alpha, beta scalars
//        (   x   )   (   0  )     x  (n-1)-dimensional vector
//
// Because of H'* H = I it follows that
//
//    alpha^2 + x^T x = beta^2    ==> beta = sqrt(alpha^2 + x^T x)
//
//    H'( beta )  = ( alpha )     ==> tau = (beta - alpha)/beta
//      (  0   )    (  x    )           v = xscal*x,
//                                  xscal = 1/(beta - alpha)
//
//    If x = 0, then tau = 0 and H = I, otherwise  1 <= tau <= 2.
//
// This is an edited version of the LAPACK routine DLARFG
//
// Variables in the calling sequence:
// ----------------------------------
// N      I   IN   Dimension of H
// ALPHA  D   IN   The scalar alpha
//            OUT  The scalar beta
// X      D   IN   The given vector x of dimension n - 1
//            OUT  The vector v
// INCX   I   IN   The increment between elements of X, INCX .NE. 0
// TAU    D   OUT  The scalar tau
// SAFMIN D   IN   Safe minimum such that 1.0/SAFMIN does not
//                 overflow
int housg(int n, double *alpha, double *x, int incx, double *tau, double safmin)
{
  const double one = 1.0, zer = 0.0;
  int ix = 0, knt = 0, kx = 0, nm1 = 0;
  double a1 = 0.0, a2 = 0.0, beta = 0.0, rsafmn = 0.0, tmp = 0.0, xnorm = 0.0,
         xscal = 0.0;

  // H is defined as the identity whenever incx==0; implies don't iterate
  // through x
  if (incx == 0)
  {
    *tau = zer;
    return -2;
  }

  // H is defined as the identity whenever n<1
  if (n < 1)
  {
    *tau = zer;
    return -1;
  }

  // H is defined as the identity whenever n=1
  if (n == 1)
  {
    *tau = zer;
    return 0;
  }

  // Dimension of x
  nm1 = n - 1;
  // Norm of x
  xnorm = dnrm2(nm1, x, incx);

  // H is the identity whenever xnorm=0 with alpha=beta and tau=0
  if (xnorm == zer)
  {
    *tau = zer;
    return 0;
  }

  // General case
  kx = 1;
  if (incx < 0) kx = 1 - (n - 1) * incx;

  // Compute alpha, beta, tau, and xscal
  a1 = xnorm;
  a2 = fabs(*alpha);
  if (a1 < a2)
  {
    a1 = a2;
    a2 = xnorm;
  }
  if (a2 == zer)
  {
    beta = a1;
  }
  else
  {
    tmp = a1 / a2;
    beta = a2 * sqrt(one + tmp * tmp);
  }
  if (*alpha > zer) beta = -beta;

  // Test for loss of accuracy
  if (fabs(beta) >= safmin)
  {
    *tau = (beta - *alpha) / beta;
    xscal = one / (*alpha - beta);
    *alpha = beta;
  }
  else  // xnorm, beta may be inaccurate; scale x and recompute
  {
    rsafmn = one / safmin;
    knt = 0;
    do
    {
      knt = knt + 1;
      if (incx == 1)
      {
        for (int i = 1; i <= nm1; i++)
        {
          x(i) = rsafmn * x(i);
        }
      }
      else
      {
        ix = kx;
        for (int i = 1; i <= nm1; i++)
        {
          x(ix) = rsafmn * x(ix);
          ix = ix + incx;
        }
      }
      beta = beta * rsafmn;
      *alpha = (*alpha) * rsafmn;
    } while (fabs(beta) < safmin);

    // new beta satisfies safmin <= beta <= 1.0
    xnorm = dnrm2(nm1, x, incx);
    a1 = xnorm;
    a2 = fabs(*alpha);
    if (a1 < a2)
    {
      a1 = a2;
      a2 = xnorm;
    }
    if (a2 == zer)
    {
      beta = a1;
    }
    else
    {
      tmp = a1 / a2;
      beta = a2 * sqrt(one + tmp * tmp);
    }
    if (*alpha > zer) beta = -beta;
    *tau = (beta - *alpha) / beta;
    xscal = one / (*alpha - beta);
    *alpha = beta;
    for (int j = 1; j <= knt; j++)
    {
      *alpha = (*alpha) * safmin;
    }
  }
  if (incx == 1)
  {
    for (int i = 1; i <= nm1; i++)
    {
      x(i) = xscal * x(i);
    }
  }
  else
  {
    ix = kx;
    for (int i = 1; i <= nm1; i++)
    {
      x(ix) = xscal * x(ix);
      ix = ix + incx;
    }
  }
  return 0;
}

// Multiplies a given M x N matrix A from the (L)eft by an
// M-dimensional Householder reflector
//
//    H = I - tau*x*x',   H' * H = I,  x = ( 1 )
//                                         ( v )
//
// that is, overwrite A by the product H * A
//
// Variables in the calling sequence:
// ----------------------------------
// M     I   IN   The number of rows of A
// N     I   IN   The number of columns of A
// X     D   IN   The given vector x. x(1) = 1.0 is enforced
// INCX  I   IN   The increment between elements of X, INCX .NE. 0
// TAU   D   IN   The scalar tau
// A     D   IN   The given matrix of dimension M x N
//           OUT  The computed matrix product H * A
// LDA   I   IN   The leading dimension of the array A, LDA >= M
// IER   I   OUT  Error indicator
//                IER = 0  no error
//                IER = 1  input-data error
int housl(int m, int n, double *x, int incx, double tau, double *a, int lda,
          int *ier)
{
  const double zer = 0.0;
  int ix = 0, kx = 0;
  double sum = 0.0, tmp = 0.0;

  *ier = 0;

  // Valid data
  if ((m < 0) || (n < 0) || (incx == 0) || (lda < m))
  {
    *ier = 1;
    return 1;
  }

  // Quick return
  if ((m == 0) || (n == 0) || (tau == 0.0))
  {
    *ier = 0;
    return 0;
  }

  kx = 1;
  if (incx < 0) kx = 1 - (m - 1) * incx;

  // Compute w = A' * x and  A := A - tau * x * w'

  for (int j = 1; j <= n; j++)
  {
    sum = a(1, j);
    if (m > 1)
    {
      ix = kx;
      for (int i = 2; i <= m; i++)
      {
        ix = ix + incx;
        sum = sum + a(i, j) * x(ix);
      }
    }
    if (sum != zer)
    {
      tmp = -tau * sum;
      a(1, j) = a(1, j) + tmp;
      if (m > 1)
      {
        ix = kx;
        for (int i = 2; i <= m; i++)
        {
          ix = ix + incx;
          a(i, j) = a(i, j) + x(ix) * tmp;
        }
      }
    }
  }
  return *ier;
}
