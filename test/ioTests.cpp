#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "io.hpp"

#define DOUBLE_NEAR(x) DoubleNear((x), 4.0 * DBL_EPSILON)

using namespace testing;
using namespace std;

class TestThatIO : public Test {};

TEST_F(TestThatIO, FindsROCSubdirectoryInCurrentWorkingDirectory)
{
  string roc{"roc"};
  ASSERT_NO_THROW(cwd_path_to(roc));
}

TEST_F(TestThatIO, ThrowsExceptionWhenItDoesNotFindGivenSubdirectoryInCurrentWorkingDirectory)
{
  string roc{"foo"};
  ASSERT_THROW(cwd_path_to(roc), pathException);
}

TEST_F(TestThatIO, ReadsData)
{
  filesystem::path dn;
  try
  {
    string roc{"roc"};
    dn = cwd_path_to(roc);
  }
  catch (exception& e)
  {
    cout << e.what() << '\n';
  }
  
  string filename{dn / "test" / "coeff.txt"};
  // cout << filename << "\n";
  EXPECT_THAT(filename, EndsWith(string("roc/test/coeff.txt")));
  
  ifstream input{filename};
  if (!input)
  {
      perror("Error opening input file");
      //return -1;
  }
  
  vector<double> coeffs;
  double time{0};
  double scale{0};
  string s;
  cout.precision(16);
  cout << scientific;
  while (getline(input, s))
  {
    coeffs = read_numbers(s);

    time = coeffs.back();
    coeffs.pop_back();

    scale = coeffs.back();
    coeffs.pop_back();

    //cout << coeffs;
    //cout << "T     = " << time << '\n';
    //cout << "Scale = " << scale << '\n';
    //cout << "Size  = " << coeffs.size() << '\n';
    EXPECT_THAT(time, Ge(0.0));
    EXPECT_THAT(scale, Gt(0.0));
    EXPECT_THAT(coeffs.size(), Eq(31));
  }
  input.close();    
}
