#ifndef sixterm_hpp
#define sixterm_hpp

#include <vector>
#include "matrix.hpp"
#include "vectorf.hpp"

double sixterm(const std::vector<double> &coeff, const double &scale,
               double &rc, double &order);
void constructSixTermSystem(const vectorf<double> &, const int,
                            matrix<double> &, vectorf<double> &);
void testBeta4(double beta4);
void testRCSix(double);
void testCosTheta(double);
double testSingularityOrder(double singularityOrder1, double singularityOrder2);

#endif
