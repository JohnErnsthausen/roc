#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>
#include <iostream>
#include <vector>

#include "matrix.hpp"

#define DOUBLE_NEAR(x) DoubleNear((x), DBL_EPSILON)

using namespace testing;
using namespace std;

class TestThatMatrix : public Test
{
 public:
  int m{3}, n{3};
  double epsilon{DBL_EPSILON};
  vector<double> a;
  vector<double> x;
  vector<double> b;

  void SetUp() override
  {
    a.assign(m * n, 0.0);
    x.assign(n, 0.0);
    b.assign(m, 0.0);
  }

  void TearDown() override {}
};

TEST_F(TestThatMatrix, CanConstructAnUninitializedMatrixSepcifyingRowsCols)
{
  ASSERT_NO_THROW(matrix<double> M(2,2));
}

TEST_F(TestThatMatrix, ConstructorRowsColsThrowsWheneverRowsIsNegative)
{
  // Interesting. Didn't know this was automatically taken care of by STL
  ASSERT_THROW(matrix<double> M(-1,2), std::bad_array_new_length);
}

TEST_F(TestThatMatrix, ConstructorRowsColsThrowsWheneverRowsIsZero)
{
  ASSERT_ANY_THROW(matrix<double> M(0,2));
}

TEST_F(TestThatMatrix, ConstructorRowsColsThrowsWheneverColsIsNegative)
{
  // Interesting. Didn't know this was automatically taken care of by STL
  ASSERT_THROW(matrix<double> M(2,-1), std::bad_array_new_length);
}

TEST_F(TestThatMatrix, ConstructorRowsColsThrowsWheneverColsIsZero)
{
  ASSERT_THROW(matrix<double> M(2,0), std::exception);
}

TEST_F(TestThatMatrix, get_rowsWorks)
{
  matrix<double> M(4, 7);

  // Validate these getters
  EXPECT_THAT(M.get_rows(), Eq(4));
}

TEST_F(TestThatMatrix, get_colsWorks)
{
  matrix<double> M(4, 7);

  // Validate these getters
  EXPECT_THAT(M.get_rows(), Eq(4));
}

TEST_F(TestThatMatrix, CanConstructAnUninitializedMatrixSepcifyingRowsColsANDInitialiseIt)
{
  matrix<int> M(4, 7);
  
  // Must be able to modify the element to do this
  M(0,0) = 1;M(0,1) = 5;M(0,2) = 9 ;M(0,3) = 13;M(0,4) = 17;M(0,5) = 21;M(0,6) = 25;  
  M(1,0) = 2;M(1,1) = 6;M(1,2) = 10;M(1,3) = 14;M(1,4) = 18;M(1,5) = 22;M(1,6) = 26; 
  M(2,0) = 3;M(2,1) = 7;M(2,2) = 11;M(2,3) = 15;M(2,4) = 19;M(2,5) = 23;M(2,6) = 27; 
  M(3,0) = 4;M(3,1) = 8;M(3,2) = 12;M(3,3) = 16;M(3,4) = 20;M(3,5) = 24;M(3,6) = 28;

  EXPECT_THAT(M.get_rows(), Eq(4));
  EXPECT_THAT(M.get_cols(), Eq(7));

  // Must be able to access, but not modify, the element to do this
  for(size_t j{0}; j < M.get_cols(); j++)
    for(size_t i{0}; i < M.get_rows(); i++)
      EXPECT_THAT(M(i,j), Eq(M.get_rows()*j + i + 1) );

  // Printing verification does work
  // cout << "M =\n"
  //      << M << '\n';
}

TEST_F(TestThatMatrix, ImplementsCopyConstructorCorrectly)
{
  matrix<double> M1(4, 2);

  M1(0,0) = 1;M1(0,1) = 5;
  M1(1,0) = 2;M1(1,1) = 6;
  M1(2,0) = 3;M1(2,1) = 7;
  M1(3,0) = 4;M1(3,1) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test the copy constructor
  matrix<double> M2{M1};

  // Test that matrix M1 hasn't changed
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test that matrix M2 equals M1 (Copy performed correctly)
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );

  // Test that copies are independent, the default copy constructor overridden
  M2(0, 1) = 3;
  
  // Test that matrix M1 hasn't changed
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test that matrix M2 has 3 in correct position
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      if( i==0 and j==1 )
      {
        EXPECT_THAT(M2(i,j), Eq(3) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );
      }
}

TEST_F(TestThatMatrix, ImplementsCopyAssignmentCorrectly)
{
  matrix<double> M1(4, 2);

  M1(0,0) = 1;M1(0,1) = 5;
  M1(1,0) = 2;M1(1,1) = 6;
  M1(2,0) = 3;M1(2,1) = 7;
  M1(3,0) = 4;M1(3,1) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test the copy constructor
  matrix<double> M2(4,2);

  // If quick return doesn't work, them no code coverage! No Idea how else to test this.
  M1 = M1;

  M2 = M1;

  // Test that matrix M1 hasn't changed
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test that matrix M2 equals M1 (Copy performed correctly)
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );

  // Test that copies are independent, the default copy constructor overridden
  M2(1, 0) = 4;
  
  // Test that matrix M1 hasn't changed
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test that matrix M2 has 3 in correct position
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      if( i==1 and j==0 )
      {
        EXPECT_THAT(M2(i,j), Eq(4) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );
      }
}

TEST_F(TestThatMatrix, ImplementsMoveConstructorCorrectly)
{
  matrix<double> M1(4, 2);

  M1(0,0) = 1;M1(0,1) = 5;
  M1(1,0) = 2;M1(1,1) = 6;
  M1(2,0) = 3;M1(2,1) = 7;
  M1(3,0) = 4;M1(3,1) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test the Move constructor. Move is from smart pointer?
  matrix<double> M2 = move(M1);

  // Test that matrix M1 hasn't changed
  EXPECT_THAT(M1.get_cols(), Eq(0));
  EXPECT_THAT(M1.get_rows(), Eq(0));
  EXPECT_EQ(&M1(0,0), nullptr);

  // Test that matrix M2 equals old M1 (Move performed correctly)
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );

  // Test that copies are independent, the default copy constructor overridden
  M2(0, 0) = 5;
  
  // Test that matrix M2 has 5 in correct position
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      if( i==0 and j==0 )
      {
        EXPECT_THAT(M2(i,j), Eq(5) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );
      }
}

TEST_F(TestThatMatrix, ImplementsMoveAssignmentCorrectly)
{
  matrix<double> M1(4, 2);

  M1(0,0) = 1;M1(0,1) = 5;
  M1(1,0) = 2;M1(1,1) = 6;
  M1(2,0) = 3;M1(2,1) = 7;
  M1(3,0) = 4;M1(3,1) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  // Test the Move constructor. Move is from smart pointer?
  matrix<double> M2(4, 2);

  // If quick return doesn't work, them matrix erases itself and this test cannot pass.
  M1 = move(M1);

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{0}; j < M1.get_cols(); j++)
    for(size_t i{0}; i < M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*j + i + 1) );

  M2 = move(M1);

  // Test that matrix M1 hasn't changed
  EXPECT_THAT(M1.get_cols(), Eq(0));
  EXPECT_THAT(M1.get_rows(), Eq(0));
  EXPECT_EQ(&M1(0,0), nullptr);

  // Test that matrix M2 equals old M1 (Move performed correctly)
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );

  // Test that copies are independent, the default copy constructor overridden
  M2(0, 0) = 5;
  
  // Test that matrix M2 has 5 in correct position
  for(size_t j{0}; j < M2.get_cols(); j++)
    for(size_t i{0}; i < M2.get_rows(); i++)
      if( i==0 and j==0 )
      {
        EXPECT_THAT(M2(i,j), Eq(5) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*j + i + 1) );
      }
}

TEST_F(TestThatMatrix, CannotDoCopyAssignmentWheneverDimensionsDonotAgree)
{
  matrix<double> M1(4, 2);
  matrix<double> M2(2, 4);

  ASSERT_THROW(M1 = M2, std::exception);
}

TEST_F(TestThatMatrix, CannotDoMoveAssignmentWheneverDimensionsDonotAgree)
{
  matrix<double> M1(4, 2);
  matrix<double> M2(2, 4);

  ASSERT_THROW(M1 = move(M2), std::exception);
}


// TEST_F(TestThatMatrix, TWO)
// {
//     //try
//     //{
//         // Constructor with two integers: create an UNINITIALIZED 3x4 matrix
//         matrix<int> A(3, 4);
//         cout << "A:\n"
//              << A;
//         // Constructor with one vector: create a 3x3 matrix with 1, 2, 3 on the diagonal
//         matrix<int> B(vector<int>{1, 2, 3});
//         cout << "B:\n"
//              << B;
//         // Constructor with one initializer_list: create a 4x4 matrix with 1, 2, 3, 4 on the diagonal
//         matrix<int> C{1, 2, 3, 4};
//         cout << "C:\n"
//              << C;
//         // Constructor with two integers and one vector: create a 2x3 matrix with the given elements
//         matrix<int> D(2, 3, vector<int>{1, 2, 3, 4, 5, 6});
//         cout << "D:\n"
//              << D;
//         // Constructor with two integers and one initializer_list: create a 2x2 matrix with the given elements
//         matrix<int> E(2, 2, {1, 2, 3, 4});
//         cout << "E:\n"
//              << E;
// 
//         // Demonstration of some of the overloaded operators
//         D(0, 2) = 7;
//         cout << "D after D(0, 2) = 7:\n"
//              << D;
//         matrix<int> F = D * B;
//         cout << "F = D * B:\n"
//              << F;
//         cout << "D + F:\n"
//              << D + F;
//         cout << "7 * B:\n"
//              << 7 * B;
// 
//         // initializer_list constructor will be used: create a 2x2 diagonal matrix with 1, 2 on the diagonal
//         cout << "matrix{1, 2}:\n";
//         cout << matrix<int>{1, 2};
//         // (size_t, size_t) constructor will be used: create an UNINITIALIZED 1x2 matrix
//         cout << "matrix(1, 2):\n";
//         cout << matrix<int>(1, 2);
//     //}
//     //catch (const matrix::index_out_of_range &e)
//     //{
//     //    cout << "Error: Matrix index out of range!\n";
//     //}
//     //catch (const matrix::zero_size &e)
//     //{
//     //    cout << "Error: Cannot create a matrix with zero rows or columns!\n";
//     //}
//     //catch (const matrix::initializer_wrong_size &e)
//     //{
//     //    cout << "Error: Initializer elements do not match the expected number of elements!\n";
//     //}
//     //catch (const matrix::incompatible_sizes_add &e)
//     //{
//     //    cout << "Error: Two matrices can only be added or subtracted if they are of the same size!\n";
//     //}
//     //catch (const matrix::incompatible_sizes_multiply &e)
//     //{
//     //    cout << "Error: Two matrices can only be multiplied if the number of columns in the first matrix is equal to the number of rows in the second matrix!\n";
//     //};  
// }
