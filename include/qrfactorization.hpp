#ifndef qrfactorization_hpp
#define qrfactorization_hpp

#include <vector>
#include "matrix.hpp"
#include "vectorf.hpp"

int qr(const int, const int, std::vector<double> &, std::vector<double> &,
       std::vector<double> &);
int qr(const int, const int, matrix<double> &, vectorf<double> &,
       vectorf<double> &);

#endif
