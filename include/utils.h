#ifndef utils_h
#define utils_h

#include <float.h>
#include <stdlib.h>
#include <tgmath.h>
#include <assert.h>

// uncomment to disable assert()
// #define NDEBUG

int sgn(double val);
double sign(double a, double b);
void swap(double *x, double *y);
int housg(int n, double *alpha, double *x, int incx, double *tau, double safemin);
double ddist2(int n, double *x, int incx, double *y, int incy, int kvec);

#endif
