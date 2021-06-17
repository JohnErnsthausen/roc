#ifndef roc_hpp
#define roc_hpp

#include <vector>

int roc(const std::vector<double> &, const double &, double &, double &);
bool convergence(double);
int derivatives2coefficients(std::vector<double> &);

#endif
