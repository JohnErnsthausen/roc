#ifndef qrfactorization_h
#define qrfactorization_h

int housg(int n, double *alpha, double *x, int incx, double *tau,
          double safemin);
int housl(int m, int n, double *x, int incx, double tau, double *a, int lda,
          int *ier);
int qrf(int m, int n, double *a, int lda, int *ipiv, double *tau, double *wrk,
        double safmin, int *ier);
int qrs(int m, int n, double *a, int lda, double *tau, double *y, double *x,
        int *ier);

#endif
