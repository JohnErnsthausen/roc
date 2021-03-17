#ifndef io_hpp
#define io_hpp

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

std::ostream &operator<<(std::ostream &out, const std::vector<double> &v);
std::vector<double> read_numbers(const std::string &in);
std::filesystem::path cwd_path_to(const std::string &name);

#endif
