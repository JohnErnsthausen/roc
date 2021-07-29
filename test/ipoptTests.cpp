#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cassert>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "Ipopt_Interface.h"
#include "nan.h"

/* A structure to store the sparsity structure of the Jacobian and Hessian
   (this is an example on how userData might be used, there is no need to
    store your matrices using this structure) */
struct _MatrixStructure {
  int nnz;
  int nzh;
  int *con_row_index;
  int *con_col_index;
  int *lhess_col_index;
  int *lhess_row_index;
};
typedef struct _MatrixStructure *MatrixStructure;

/* Constructor for MatrixStructure */
MatrixStructure MatrixStructure_Create(int nnz, int nzh)
{
  MatrixStructure MS;

  MS = (MatrixStructure)malloc(sizeof(struct _MatrixStructure));
  if (MS == NULL) {
    printf("Could not allocate the IpoptJava structure. Exiting...\n");
    exit(2);
  }

  MS->nnz = nnz;
  MS->nzh = nzh;

  if (nnz>0) {
    MS->con_row_index = (int*)malloc(nnz*sizeof(int));
    if (MS->con_row_index == NULL) {
      printf("Could not allocate space for con_row_index. Exiting...\n");
      exit(2);
    }
    MS->con_col_index = (int*)malloc(nnz*sizeof(int));
    if (MS->con_col_index == NULL) {
      printf("Could not allocate space for con_col_index. Exiting...\n");
      exit(2);
    }
  }
  if (nzh>0) {
    MS->lhess_col_index = (int*)malloc(nzh*sizeof(int));
    if (MS->lhess_col_index == NULL) {
      printf("Could not allocate space for lhess_col_index. Exiting...\n");
      exit(2);
    }
    MS->lhess_row_index = (int*)malloc(nzh*sizeof(int));
    if (MS->lhess_row_index == NULL) {
      printf("Could not allocate space for lhess_row_index. Exiting...\n");
      exit(2);
      }
  }

  return MS;
}

/* Destructor for MatrixStructure */
void MatrixStructure_Destroy(MatrixStructure MS)
{
  if(MS->nnz>0) {
    free(MS->con_row_index);
    free(MS->con_col_index);
  }
  if (MS->nzh>0) {
    free(MS->lhess_col_index);
    free(MS->lhess_row_index);
  }
  free(MS);
}

/*
  Evaluate objective function

  n:   number of optimization variables
  x:   values of optimization variables at which objective function
       is to be evaluated
  f:   value of the objective function
*/
int Eval_F(int n, double *x, double *f, void *userData)
{
  assert(x && f && userData);
  init_with_nans(1, f);

  double *X = x;
  X--;   /* This is just a trick to start counting the variables with 1
	    and not 0 */
  *f =  X[1] * X[4] * (X[1] + X[2] + X[3]) + X[3];


  int check = check_for_nan(1, f);
  if (check >= 0)
    throw std::logic_error("IPOPT objective evaluation is NaN\n");
  return 0;
}

/*
  Evaluate gradient of objective function

  n:   number of optimization variables
  x:   values of optimization variables at which gradient
       is to be evaluated
  g:   values of the gradient
*/
int Eval_G(int n, double *x, double *g, void *userData)
{
  assert(x && g && userData);
  init_with_nans(n, g);

  double *G = g;
  double *X = x;
  G--;
  X--;

  G[1] = X[4] * (2e0 * X[1] + X[2] + X[3]);
  G[2] = X[1] * X[4];
  G[3] = X[1] * X[4] + 1e0;
  G[4] = X[1] * (X[1] + X[2] + X[3]);
  G[5] = 0e0;

  int check = check_for_nan(n, g);
  if (check >= 0)
    throw std::logic_error("IPOPT gradient of objective evaluation has NaN(s)\n");
  return 0;
}

/*
  Evaluate constraints

  n:   number of optimization variables
  x:   values of optimization variables at which constraints
       are to be evaluated
  m:   number of constraints
  c:   values of the constraints
*/
int Eval_C(int n, double *x, int m, double *c, void *userData)
{
  assert(x && c && userData);
  init_with_nans(m, c);

  double *C = c;
  double *X = x;
  C--;
  X--;

  C[1] = X[1] * X[2] * X[3] * X[4] - X[5] - 25e0;
  C[2] = X[1] * X[1] + X[2] * X[2] + X[3] * X[3] + X[4] * X[4] - 40e0;

  int check = check_for_nan(m, c);
  if (check >= 0)
    throw std::logic_error("IPOPT constraint evaluation has NaN(s)\n");
  
  return 0;
}

/*
  Evaluate constraint Jacobian

  task: if 0, then only return number of nonzeros in nz
        otherwise:

  n:    number of optimization variables
  x:    values of optimization variables at which constraint Jacobian
        is to be evaluated
  nz:   number of nonzeros in Jacobian
  A:    values of the Jacobian entries
  Arow: row indices of the Jacobian entries
  Acol: column indices of the Jacobian entries
*/
int Eval_A(int task,
           int n,
           double *x,
           int *nz,
           double *A,
           int *Arow,
           int *Acol,
           void *userData)
{
  assert(x && nz && A && Arow && Acol && userData);
  //if (task != 0) init_with_nans(*nz, A);

  double *A_ = A;
  double *X = x;
  int i;
  MatrixStructure MS = (MatrixStructure)userData;

  A_--;
  X--;

  if (task == 0) {
    *nz = 9;   /* assign number of nonzeros */
  }
  else {
    for (i = 0; i < *nz; i++) {
      Arow[i] = MS->con_row_index[i]; /* rows */
      Acol[i] = MS->con_col_index[i]; /* cols */
    }
    A_[1] = X[2]        * X[3] * X[4];
    A_[2] = X[1]        * X[3] * X[4];
    A_[3] = X[1]        * X[2] * X[4];
    A_[4] = X[1]        * X[2] * X[3];
    A_[5] = -1e0;
    A_[6] = 2e0 * X[1];
    A_[7] = 2e0 * X[2];
    A_[8] = 2e0 * X[3];
    A_[9] = 2e0 * X[4];
  }

  //int check = check_for_nan(*nz, A);
  //if (check >= 0)
  //  throw std::logic_error("IPOPT constraint Jacobian evaluation has NaN(s)\n");
  
  return 0;
}

/*
  Evaluate the Hessian of the Lagrangian function

  task:   if 0, the only return number of nonzeros in nnzh
          otherwise:

  n:      number of optimization variables
  x:      values of optimization variables at which constraint Jacobian
          is to be evaluated
  m:      number of constraints
  lambda: values of multipliers
  nnzh:   number of nonzeros in Jacobian
  hess:   values of the Hessian entries
  irnh:   row indices of the Hessian entries
  icnh:   column indices of the Hessian entries
*/
int Eval_H(int task, int n, double *x, int m,
           double *lambda, int *nnzh, double *hess,
           int *irnh, int *icnh, void *userData)
{
  double *LAM = lambda;
  double *X = x;
  double *HESS = hess;
  int i;
  MatrixStructure MS = (MatrixStructure)userData;

  LAM--;
  X--;
  HESS--;

  if (task == 0) {
    *nnzh = 10;
  }
  else {
    *nnzh = 10;
    for (i = 0; i < *nnzh; i++) {
      irnh[i] = MS->lhess_row_index[i];
      icnh[i] = MS->lhess_col_index[i];
      hess[i] = 0.0; /* this is critical here, because the lines below are
			only adding things */
    }

    /*
     * Objective function
     */
    HESS[1] = 2e0 * X[4];
    HESS[2] = X[4];
    HESS[4] = X[4];
    HESS[7] = 2e0 * X[1] + X[2] + X[3];
    HESS[8] = X[1];
    HESS[9] = X[1];

    /*
     * first constraint
     */
    HESS[2] += LAM[1] * X[3] * X[4];
    HESS[4] += LAM[1] * X[2] * X[4];
    HESS[5] += LAM[1] * X[1] * X[4];
    HESS[7] += LAM[1] * X[2] * X[3];
    HESS[8] += LAM[1] * X[1] * X[3];
    HESS[9] += LAM[1] * X[1] * X[2];

    /*
     * second constraint
     */ 
    HESS[1] += LAM[2] * 2e0;
    HESS[3] += LAM[2] * 2e0;
    HESS[6] += LAM[2] * 2e0;
    HESS[10]+= LAM[2] * 2e0;

  }

  return 0;
}

using namespace testing;

class TestThatIPOPTSixTermAnalysis : public Test
{
 public:
  double epsilon{DBL_EPSILON};

  void SetUp() override {}

  void TearDown() override {}
};

TEST_F(TestThatIPOPTSixTermAnalysis, works)
{
  /* Those first variables define the problem size: */
  fint n = 5;      /* number of optimization variables */
  fint m = 2;      /* number of equality constraints */
  fint nlb = 5;    /* number of lower bounds on variables */
  fint nub = 4;    /* number of upper bounds on variables */

  int i;

  /* The indices for the lower bounds */
  fint ilb[] = {
    1,   2,   3,   4,   5
  };

  /* The values for the lower bounds */
  real bnds_l[] = {
    1e0, 1e0, 1e0, 1e0, 0e0
  };

  /* The indicies for the upper bounds */
  fint iub[] = {
    1,   2,   3,   4
  };

  /* The values for the upper bounds */
  real bnds_u[] = {
    5e0, 5e0, 5e0, 5e0
  };

  /* Starting point */
  real xstart[] = {
    1e0, 5e0, 5e0, 1e0, -24e0
  };

  /* number of nonzeros in the constraint Jacobian */
  fint nnz = 9;

  /* row indices of the constraint Jacobian */
  fint con_row_index[] = { /* acon */
    1, 1, 1, 1, 1, 2, 2, 2, 2
  };

  /* column indices of the constraint Jacobian */
  fint con_col_index[] = { /* avar */
    1, 2, 3, 4, 5, 1, 2, 3, 4
  };

  /* number of nonzeros in the (symmetric!) Largrangian Hessian */
  fint nzh = 10;

  /* row indices for the lower diagonal part of the Largrangian Hessian */
  fint lhess_row_index[] = { /* irnh */
    1, 2, 2, 3, 3, 3, 4, 4, 4, 4
  };

  /* column indices for the lower diagonal part of the Largrangian Hessian */
  fint lhess_col_index[] = { /* icnh */
    1, 1, 2, 1, 2, 3, 1, 2, 3, 4
  };

  real Obj;
  real *x=NULL;
  real *v_lb=NULL;
  real *v_ub=NULL;
  real *lambda=NULL;
  real *c=NULL;
  fint iter, retval;

  MatrixStructure MS;

  Ipopt problem = NULL;

  /* First create an Ipopt object, which will store certain quantities,
     such as size of the problem, values of the variables (initial point, and
     later the solution), etc
  */
  problem = Ipopt_Create(n, m, nlb, ilb, bnds_l, nub, iub, bnds_u, Eval_F,
			 Eval_C, Eval_G, Eval_A, Eval_H);

  if (problem == NULL) {
    printf("Could not allocate Ipopt problem\n");
    exit(2);
  }

  MS = MatrixStructure_Create(nnz, nzh);

  /* Take care of constraint Jacobian structure (used in Eval_A) */
  for (i = 0; i < nnz; i++ ) {
    MS->con_col_index[i] = con_col_index[i];
    MS->con_row_index[i] = con_row_index[i];
  }

  /* Take care of Lagrangian Hessian structure (used in Eval_H) */
  for (i = 0; i < nzh; i++ ) {
    MS->lhess_col_index[i] = lhess_col_index[i];
    MS->lhess_row_index[i] = lhess_row_index[i];
  }

  /*
   * Set the algorithmic parameters (see description of INIT_PARAMS)
   */
  char dtol_s[] = "DTOL", iprint_s[] = "IPRINT", ioutput_s[] = "IOUTPUT",
       ifile_s[] = "IFILE", ifull_s[] = "IFULL", iquasi_s[] = "IQUASI",
       isymsolv_s[] = "ISYMSOLV";
  Ipopt_AddParam(problem, dtol_s, 1.0e-8);
  Ipopt_AddParam(problem, iprint_s, -1);
  Ipopt_AddParam(problem, ioutput_s, 1);
  Ipopt_AddParam(problem, ifile_s, 0);
  Ipopt_AddParam(problem, ifull_s, 1);
  Ipopt_AddParam(problem, iquasi_s, 6);
  Ipopt_AddParam(problem, isymsolv_s, 1);

  /*
    Get the memory for the result from IPOPT
  */
  x = (real*)malloc(sizeof(real)*n);
  if( nlb>0 ) v_lb = (real*)malloc(sizeof(real)*nlb);
  if( nub>0 ) v_ub = (real*)malloc(sizeof(real)*nub);
  if( m>0 ) {
    c = (real*)malloc(sizeof(real)*m);
    lambda = (real*)malloc(sizeof(real)*m);
  }

  /*
    Set the starting point
  */
  for( i=0; i<n; i++) x[i] = xstart[i];

  /*
   * Call the optimization routine
   */
  retval = Ipopt_Solve(problem, (void*)MS, x, &Obj, v_lb, v_ub, lambda,
		       c, &iter);

  /*
   * Output:
   */
  if (retval == 0) {
    printf("\nThe final values of the objective function: %15.8e\n",Obj);

    printf("\nThe optimal values of X are:\n\n");
    for (i = 0; i < n; i++) {
      printf("X  (%5d) = %15.8e\n", i+1, x[i]);
    }

    printf("\nThe multipliers for the lower bounds are:\n\n");
    for (i = 0; i < nlb; i++) {
      printf("V_L(%5d) = %15.8e\n", ilb[i], v_lb[i]);
    }

    printf("\nThe multipliers for the upper bounds are:\n\n");
    for (i = 0; i < nub; i++) {
      printf("V_U(%5d) = %15.8e\n", iub[i], v_ub[i]);
    }

    printf("\nThe multipliers for the equality constraints are:\n\n");
    for (i = 0; i < m; i++) {
      printf("LAM(%5d) = %15.8e\n", i+1, lambda[i]);
    }
    printf("\n");
  }
  else {
    printf("\nAn error has occurred after %d iterations\n", iter);
    printf("The error code is %d\n\n", retval);
    exit(2);
  }

  /*
   * Clean up and return
   */
  if( m>0 ) {
    free(lambda);
    free(c);
  }
  if( nub>0 ) free(v_ub);
  if( nlb>0 ) free(v_lb);
  free(x);
  Ipopt_Destroy(problem);
  MatrixStructure_Destroy(MS);
}

