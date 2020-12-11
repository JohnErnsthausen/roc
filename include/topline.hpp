#ifndef topline_hpp
#define topline_hpp

#include <vector>
#include "matrix.hpp"
#include "vectorf.hpp"

double topline(const std::vector<double> &, const double &, double &, double &);
void constructLinearLeastSquaresSystem(const std::vector<double> &, const int,
                                       matrix<double> &, vectorf<double> &);

#endif
