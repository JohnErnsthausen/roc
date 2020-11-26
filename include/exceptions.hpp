#ifndef exceptions_hpp
#define exceptions_hpp

#include <exception>
#include <string>

class sayMessage : public std::exception
{
 public:
  explicit sayMessage(const std::string &message) : message_(message){};
  const char *what() const noexcept override { return message_.c_str(); }

 private:
  std::string message_;
};

#endif
