#ifndef topline_hpp
#define topline_hpp

#include <vector>
#include "matrix.hpp"
#include "vectorf.hpp"

void constructLinearLeastSquaresRow(const std::vector<double> &, const int,
                                    double &, double &, double &);
void constructLinearLeastSquaresSystem(const std::vector<double> &, const int,
                                       matrix<double> &, vectorf<double> &);
double errorTopLine(const std::vector<double> &, const vectorf<double> &);
double topline(const std::vector<double> &, const double &, double &, double &);

#endif
