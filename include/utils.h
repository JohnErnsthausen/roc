#ifndef utils_h
#define utils_h

#include <float.h>
#include <stdlib.h>
#include <tgmath.h>
#include <assert.h>

#include <stdio.h>

// uncomment to disable assert()
// #define NDEBUG

int sgn(double val);
double sign(double a, double b);
void dswap(double *x, double *y);
void iswap(int *x, int *y);

int housg(int n, double *alpha, double *x, int incx, double *tau, double safemin);
int housl(int m, int n, double *x, int incx, double tau, double *a, int lda, int *ier);
double diff_avoids_subtractive_cancellation(double dx, double dy);
int qrf( int m, int n, double *a, int lda, int *ipiv, double *tau,
         double *wrk, double safmin, int *ier );
int qrs( int m, int n, double *a, int lda, double *tau, double *y, double *x, int *ier );

double ddist2(int n, double *x, int incx, double *y, int incy, int kvec);

#endif
