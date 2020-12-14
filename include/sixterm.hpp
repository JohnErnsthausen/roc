#ifndef sixterm_hpp
#define sixterm_hpp

#include <vector>

double sixterm(const std::vector<double> &coeff, const double &scale,
               double &rc, double &order);
void testBeta4(double beta4);
void testRCSix(double);
void testCosTheta(double);
double testSingularityOrder(double singularityOrder1, double singularityOrder2);

#endif
