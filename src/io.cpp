#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "exceptions.hpp"
#include "io.hpp"

// Overloaded operator << to print out a formatted vector
std::ostream &operator<<(std::ostream &out, const std::vector<double> &v)
{
  size_t s{v.size() - 1};
  out << "(";
  for (size_t i{0}; i < s; i++) out << v[ i ] << ", ";
  out << v[ s ] << ")\n";
  return out;
}

// Function to convert a string of space-separated numbers into a vector
std::vector<double> read_numbers(const std::string &in)
{
  // A vector to store the converted numbers
  std::vector<double> v;
  // A string to store each substring we read with getline
  std::string s;
  // Create an input string stream based on the function's argument
  std::istringstream string_stream(in);
  // Read space-separated values with getline
  while (std::getline(string_stream, s, ' '))
  {
    // Convert each substring to a double with stod, and store it in the vector
    v.push_back(stod(s));
  }
  return v;
}

// Fully qualified pathname to repository expected to be a parent directory
// of the current path
std::filesystem::path cwd_path_to(const std::string &name)
{
  std::filesystem::path cwd = std::filesystem::current_path();
  std::filesystem::path dn;
  for (auto it = cwd.begin(); it != cwd.end(); ++it)
  {
    dn /= *it;
    if (*it == name) return dn;
  }
  std::string message{"Folder is not a subdirectory of the current path"};
  throw sayMessage(message);
}
