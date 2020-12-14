#ifndef threeterm_hpp
#define threeterm_hpp

#include <vector>
#include "matrix.hpp"
#include "vectorf.hpp"

double threeterm(const std::vector<double> &, const double &, double &,
                 double &);
void constructThreeTermSystem(const vectorf<double> &, const int,
                              matrix<double> &, vectorf<double> &);
void testRCThree(double);
void testOrder(double);

#endif
