#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "io.hpp"
#include "roc.hpp"

int main()
{
  // Find the input file filename in predifined and know directory structure of
  // roc
  std::filesystem::path dn;
  try
  {
    std::string roc{"roc"};
    dn = cwd_path_to(roc);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << '\n';
  }

  std::string filename{dn / "test" / "coeff.txt"};
  std::cout << filename << "\n";

  std::ifstream input{filename};
  if (!input)
  {
    perror("Error opening input file");
    return -1;
  }

  // Variables to interact with roc
  std::vector<double> coeffs;
  double time{0};
  double scale{0};
  double rc{0};
  double order{0};

  // variable for reading lines
  std::string s;

  // Output format
  std::cout.precision(16);
  std::cout << std::scientific;

  // Read the file and interact with roc
  while (getline(input, s))
  {
    coeffs = read_numbers(s);

    time = coeffs.back();
    coeffs.pop_back();

    scale = coeffs.back();
    coeffs.pop_back();

    // Print out input
    std::cout << coeffs;
    std::cout << "T     = " << time << '\n';
    std::cout << "Scale = " << scale << '\n';
    std::cout << "Size  = " << coeffs.size() << '\n';

    // Call ROC here
    roc(coeffs, scale, rc, order);
  }
  // Close the input file
  input.close();
}
