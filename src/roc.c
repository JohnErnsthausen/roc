#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "qrfactorization.h"

#define map(i, j) ((j)-1) * lda + ((i)-1)
#define a(i, j) a[ map(i, j) ]
#define map1(i) (i) - 1
#define x(i) x[ map1(i) ]
#define y(i) y[ map1(i) ]
#define ipiv(i) ipiv[ map1(i) ]
#define tc(i) tc[ map1(i) ]

//  error code 1 if less than 10 coefficients to estimate Rc

int roc(int num_tc, double *tc, double scale, int kstart, double *rc,
        double *slope, double *intercept)
{
  if (num_tc - kstart < 11)
  {
    printf(
        "Message from ROC: Need at least 10 Taylor coefficients to estimate "
        "Rc\n");
    return 1;
  }

  int m = num_tc - kstart + 1;
  int lda = m, n = 2, min_m_n = 0;
  int res = 0, ier = 0;
  double *a = NULL;
  int *ipiv = NULL;
  double *tau = NULL;
  double *wrk = NULL;
  double *x = NULL;
  double *y = NULL;
  double safmin = DBL_MIN;

  min_m_n = (m < n) ? m : n;

  a = (double *)calloc((size_t)(lda * n), (size_t)sizeof(double));
  if (a == NULL)
  {
    printf("Failed to allocate matrix A memory!\n");
    assert(a);
  }

  tau = (double *)calloc((size_t)min_m_n, (size_t)sizeof(double));
  if (tau == NULL)
  {
    printf("Failed to allocate tau memory!\n");
    assert(tau);
  }

  wrk = (double *)calloc((size_t)n, (size_t)sizeof(double));
  if (wrk == NULL)
  {
    printf("Failed to allocate wrk memory!\n");
    assert(wrk);
  }

  ipiv = (int *)calloc((size_t)n, (size_t)sizeof(int));
  if (ipiv == NULL)
  {
    printf("Failed to allocate ipiv memory!\n");
    assert(ipiv);
  }

  x = (double *)calloc((size_t)n, (size_t)sizeof(double));
  if (x == NULL)
  {
    printf("Failed to allocate space for result x!\n");
    assert(x);
  }

  y = (double *)calloc((size_t)m, (size_t)sizeof(double));
  if (y == NULL)
  {
    printf("Failed to allocate space for result y!\n");
    assert(y);
  }

  // Construct the least squares linear system
  //
  // We want to fit the tail of the Taylor series, ignore the first kstart terms
  // of the TCs.
  for (int k = kstart; k <= num_tc; k++)
  {
    int row = k - kstart + 1;
    a(row, 1) = 1.0;
    a(row, 2) = (double)k;
    y(row) = log10(fabs(tc(k)));
  }

  // Find the QRFactorization of the least squares linear system
  res = qrf(m, n, a, lda, ipiv, tau, wrk, safmin, &ier);
  if (res != 0 || ier != 0)
  {
    printf("Error QR factoring A with res[%d]!=0 or ier[%d]!=0\n", res, ier);
  }

  // Find the least squares solution of the linear system via the
  // QRFactorization
  //
  // Solve on the first (two) n rows of R as in
  // transpose(Q)*b=transpose(Q)*A*E*E*x=transpose(Q)*Q*R*E*x=R*E*x
  res = qrs(m, n, a, lda, tau, y, x, &ier);
  if (res != 0 || ier != 0)
  {
    printf(
        "Error finding the Least Squares Solution with res[%d]!=0 or "
        "ier[%d]!=0\n",
        res, ier);
  }

  // If Taylor coefficients TC(i) are unscaled, then the scaled TCs are
  // T(i) = (1/scale)^i TC(i).
  //
  // Multiply Rc by scale to get Rc when computing with scaled TCs.
  //
  // See this from ratio test as Rc is 1/L where
  //
  //   L = (1/scale)*limit of abs(TC(i+1)/TC(i)) = limit of abs(T(i+1)/T(i)).
  //
  // To justify the formula for Rc, consider the problem of fitting m and b
  //
  // log10(abs(T(i))) = y(i) = b + m*i = [1 i](b,m)
  //
  // by linear least squares. Then
  //
  // Log10(L) = Log10(abs(T(i+1))/abs(T(i))) =
  // Log10(abs(T(i+1)))-Log10(abs(T(i)))
  //          = m*(i+1)+b - (m*i+b) = m
  //
  // It follows that
  //
  // rc = 1/L = scale/pow(10, m).
  *rc = scale / pow(10, x(ipiv(2)));
  *intercept = x(ipiv(1));
  *slope = x(ipiv(2));

  free(y);
  free(x);
  free(ipiv);
  free(wrk);
  free(tau);
  free(a);

  return 0;
}
