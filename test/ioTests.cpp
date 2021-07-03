#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "io.hpp"

using namespace testing;

TEST(TestThatIO, FindsROCSubdirectoryInCurrentWorkingDirectory)
{
  std::cout << std::filesystem::current_path() << "\n";
  std::string roc{"roc"};
  ASSERT_NO_THROW(cwd_path_to(roc));
}

TEST(
    TestThatIO,
    ThrowsExceptionWhenItDoesNotFindGivenSubdirectoryInCurrentWorkingDirectory)
{
  std::string roc{"foo"};
  ASSERT_THROW(cwd_path_to(roc), std::exception);
}

TEST(TestThatIO, ReadsData)
{
  std::filesystem::path dn;
  // try
  //{
  std::string roc{"roc"};
  dn = cwd_path_to(roc);
  //}
  // catch (const std::exception& e)
  //{
  //  cout << e.what() << '\n';
  //}

  std::string filename{dn / "test" / "coeff.txt"};
  // cout << filename << "\n";
  EXPECT_THAT(filename, EndsWith(std::string("roc/test/coeff.txt")));

  std::ifstream input{filename};
  // if (!input)
  //{
  //  perror("Error opening input file");
  //  // return -1;
  //}

  std::vector<double> coeffs;
  double time{0};
  double scale{0};
  std::string s;
  std::cout.precision(16);
  std::cout << std::scientific;
  int i{0};
  while (getline(input, s))
  {
    coeffs = read_numbers(s);

    time = coeffs.back();
    coeffs.pop_back();

    scale = coeffs.back();
    coeffs.pop_back();

    // For testing exercise output operator << on first line read
    if (i == 0)
    {
      std::cout << coeffs;
      std::cout << "T     = " << time << '\n';
      std::cout << "Scale = " << scale << '\n';
      std::cout << "Size  = " << coeffs.size() << '\n';
      i++;
    }
    EXPECT_THAT(time, Ge(0.0));
    EXPECT_THAT(scale, Gt(0.0));
    EXPECT_THAT(coeffs.size(), Eq(31));
  }
  input.close();
}
