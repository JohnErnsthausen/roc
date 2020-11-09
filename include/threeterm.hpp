#ifndef threeterm_hpp
#define threeterm_hpp

#include <exception>
#include <vector>

#define TOL 1.0e-1
#define MINTERMS 10

class morecoefficients : public std::exception
{
 public:
  virtual const char *what() const throw()
  {
    return "More coefficients required. See MINTERMS.";
  }
};

class dividebyzero : public std::exception
{
 public:
  virtual const char *what() const throw() { return "Division by zero."; }
};

int threeTerm(const std::vector<double> &coeff, const double &scale, double &rc,
              double &order);

#endif
