
extern "C"
{
  // From LAPACK
  double dnrm2_(int* m, double* x, int* incx);
}

double norm2(int m, double* x)
{
  int incx{1};
  return dnrm2_(&m, x, &incx);
}

