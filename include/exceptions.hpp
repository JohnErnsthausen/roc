#ifndef exceptions_hpp
#define exceptions_hpp

#include <exception>
#include <string>

class sayMessage : public std::exception
{
 public:
  explicit sayMessage(const std::string& message) : message_(message){};
  virtual const char* what() const noexcept { return message_.c_str(); }

 private:
  std::string message_;
};

#endif
