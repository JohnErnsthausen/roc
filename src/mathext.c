#include <tgmath.h>

int sgn(double val) { return (int)(((double) 0 < val) - (val < (double) 0)); }

double sign(double a, double b) { return (b >= 0.0) ? fabs(a) : -fabs(a); }
