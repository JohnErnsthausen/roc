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

int derivatives2coefficients(std::vector<double> &coeffs)
{
  double fact=1.0;
  for(long unsigned int norder=0; norder<coeffs.size(); norder++)
  {
    coeffs[norder] = coeffs[norder]/fact;
    fact = double(norder+1)*fact;
  }
  return 0;
}

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
  int fail_count{0};
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
    std::cout << "3TA rc[" << rc << "] order [" << int(order) << "] err [" << err
              << "]\n";
  }
  catch (const std::exception &e)
  {
    fail_count++;
    std::cout << "3TA failed" << '\n';
  }
  try
  {
    err = sixterm(coeffs, scale, rc, order);
    std::cout << "6TA rc[" << rc << "] order [" << int(order) << "] err [" << err
              << "]\n";
  }
  catch (const std::exception &e)
  {
    fail_count++;
    std::cout << "6TA failed" << '\n';
  }
  try
  {
    err = topline(coeffs, scale, rc, order);
    std::cout << "TLA rc[" << rc << "] order [" << order << "] err [" << err
              << "]\n";
  }
  catch (const std::exception &e)
  {
    fail_count++;
    std::cout << "TLA failed" << '\n';
  }

  if (fail_count == 3)
  {
    std::cout << "Each method failed to resolve the RC" << '\n';
    exit(1);
  }

  return 0;
}
