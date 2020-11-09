#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "io.hpp"

using namespace std;

// Overloaded operator << to print out a formatted vector
ostream &operator<<(ostream &out, const vector<double> &v)
{
  size_t s{v.size() - 1};
  out << "(";
  for (size_t i{0}; i < s; i++)
      out << v[i] << ", ";
  out << v[s] << ")\n";
  return out;
}

// Function to convert a string of space-separated numbers into a vector
vector<double> read_numbers(const string &in)
{
  // A vector to store the converted numbers
  vector<double> v;
  // A string to store each substring we read with getline
  string s;
  // Create an input string stream based on the function's argument
  istringstream string_stream(in);
  // Read space-separated values with getline
  while (getline(string_stream, s, ' '))
  {
    // Convert each substring to a double with stod, and store it in the vector
    v.push_back(stod(s));
  }
  return v;
}

// Fully qualified pathname to repository expected to be a parent directory
// of the current path
filesystem::path cwd_path_to(const string &name)
{
  filesystem::path cwd = filesystem::current_path();
  filesystem::path dn;
  for( auto it = cwd.begin(); it != cwd.end(); ++it )
  {
    dn /= *it;
    if( *it == name ) return dn;
  }
  throw pathException();
}
