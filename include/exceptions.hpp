#ifndef exceptions_hpp
#define exceptions_hpp

#include <exception>
#include <string>

class sayMessage: public std::exception
{
public:
  explicit sayMessage(const std::string& message): message_(message) { };
  const char* what() const noexcept override { return message_.c_str(); }
private:
  std::string message_;
};

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

class QRFactorization : public std::exception
{
 public:
  virtual const char *what() const throw() { return "Error from QRFactorization"; }
};

#endif
