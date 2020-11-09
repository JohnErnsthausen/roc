#ifndef io_hpp
#define io_hpp

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

class pathException: public std::exception
{
public:
  virtual const char* what() const throw()
  {
    return "Folder is not a subdirectory of the current path";
  }
};

std::ostream &operator<<(std::ostream &out, const std::vector<double> &v);
std::vector<double> read_numbers(const std::string &in);
std::filesystem::path cwd_path_to(const std::string &name);

#endif
