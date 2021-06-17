
#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <cassert>

//#include "timings.h"

//namespace daets
//{
extern "C"
{
  // linear solver from LAPACK
  void dgels_(char* trans, int* m, int* n, int* nrhs, double* A, int* lda,
              double*, int* ldb, double* work, int* lwork, int* info);
}

void MinNormSolution(int m, int n, double* A, double* b)
{
  // Computes a least squares solutution to A*x = b. The result is
  // returned in b.
  //START_TIMING(_DGELS);
  char trans = 'N';  // no transpose
  int lda = std::max(m, 1);
  int ldb = std::max(std::max(1, n), m);
  int nrhs = 1;
  int info;

  // Determine optimal work array
  int lwork = -1;
  double work_opt;
  dgels_(&trans, &m, &n, &nrhs, A, &lda, b, &ldb, &work_opt, &lwork, &info);
  if (info < 0)
    {
      fprintf(stderr, "DGELS: the %d-th argument had an illegal value\n",
              -info);
      throw std::logic_error("Exception in DGELS");
    }

  lwork = (int)work_opt;
  
  // and allocate it
  double* pwork = new double[ lwork ];

  // Solve
  dgels_(&trans, &m, &n, &nrhs, A, &lda, b, &ldb, pwork, &lwork, &info);

  delete[] pwork;
  //END_TIMING(_DGELS);
  if (info < 0)
    {
      fprintf(stderr, "DGELS: the %d-th argument had an illegal value\n",
              -info);
      throw std::logic_error("Exception in DGELS");
    }
  if (info > 0)
    {
      fprintf(stderr, "The diagonal element %i of the triangular factor ",
              info);
      fprintf(stderr, "of A is zero, so that A does not have full rank;\n");
      fprintf(stderr, "the least squares solution could not be computed.\n");
      throw std::logic_error("Exception in DGELS");
    }
}

//}  // namespace daets
