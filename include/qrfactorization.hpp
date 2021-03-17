#ifndef qrfactorization_hpp
#define qrfactorization_hpp

#include <vector>
#include "matrix.hpp"
#include "vectorf.hpp"

int qr(const int, const int, matrix<double> &, vectorf<double> &,
       vectorf<double> &);
int factor(const int, const int, matrix<double> &, vectorf<double> &,
           vectorf<int> &);
int solve(const int, const int, matrix<double> &, vectorf<double> &,
          vectorf<double> &, vectorf<double> &);
int permute(const int, vectorf<double> &, vectorf<int> &, vectorf<double> &);

#endif
