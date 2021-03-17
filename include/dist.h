#ifndef dist_h
#define dist_h

double ddist2(int n, double *x, int incx, double *y, int incy, int kvec);
double dnrm2(int n, double *x, int incx);
void add_next_element_squared(double xi, double *sum, double *scale);

#endif
