#ifndef roc_hpp
#define roc_hpp

#include <vector>

#define MINTERMS 10
#define TOL 1e-10

int roc(const std::vector<double> &coeff, const double &scale, double &rc,
        double &order);

#endif
