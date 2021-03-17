#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "data.hpp"
#include "exceptions.hpp"
#include "io.hpp"
#include "sixterm.hpp"
#include "threeterm.hpp"
#include "topline.hpp"

bool convergence(double err)
{
  if (err > TOL)
  {
    std::cout.precision(16);
    std::cout << std::scientific;
    std::cout << "err = " << err << "\n";
  }

  return err < TOL;
}

int roc(const std::vector<double> &coeffs, const double &scale, double &rc,
        double &order)
{
  double err;

  // Check for sufficient data to perform analysis
  int k = (int)coeffs.size();
  if (k < MINTERMS)
  {
    std::string message =
        "The number of coefficients [" + std::to_string((int)coeffs.size()) +
        "] must be greater than [" + std::to_string(MINTERMS) + "].\n";
    throw sayMessage(message);
  }

  // Output format
  std::cout.precision(16);
  std::cout << std::scientific;

  // Analyse coefficients with all three analysis
  try
  {
    err = threeterm(coeffs, scale, rc, order);
    std::cout << "3TA rc[" << rc << "] order [" << order << "] err [" << err
              << "]\n";
  }
  catch (const std::exception &e)
  {
    std::cout << e.what() << '\n';
  }
  try
  {
    err = sixterm(coeffs, scale, rc, order);
    std::cout << "6TA rc[" << rc << "] order [" << order << "] err [" << err
              << "]\n";
  }
  catch (const std::exception &e)
  {
    std::cout << e.what() << '\n';
  }
  try
  {
    err = topline(coeffs, scale, rc, order);
    std::cout << "TLA rc[" << rc << "] order [" << order << "] err [" << err
              << "]\n";
  }
  catch (const std::exception &e)
  {
    std::cout << e.what() << '\n';
  }

  return 0;
}
