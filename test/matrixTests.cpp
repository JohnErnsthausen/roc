#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <initializer_list>
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
  void SetUp() override {}
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
  M(1,1) = 1;M(1,2) = 5;M(1,3) = 9 ;M(1,4) = 13;M(1,5) = 17;M(1,6) = 21;M(1,7) = 25;  
  M(2,1) = 2;M(2,2) = 6;M(2,3) = 10;M(2,4) = 14;M(2,5) = 18;M(2,6) = 22;M(2,7) = 26; 
  M(3,1) = 3;M(3,2) = 7;M(3,3) = 11;M(3,4) = 15;M(3,5) = 19;M(3,6) = 23;M(3,7) = 27; 
  M(4,1) = 4;M(4,2) = 8;M(4,3) = 12;M(4,4) = 16;M(4,5) = 20;M(4,6) = 24;M(4,7) = 28;

  EXPECT_THAT(M.get_rows(), Eq(4));
  EXPECT_THAT(M.get_cols(), Eq(7));

  // Must be able to access, but not modify, the element to do this
  for(size_t j{1}; j <= M.get_cols(); j++)
    for(size_t i{1}; i <= M.get_rows(); i++)
      EXPECT_THAT(M(i,j), Eq(M.get_rows()*(j-1) + i) );

  // // Printing verification does work
  // cout << "M =\n"
  //      << M << '\n';
}

TEST_F(TestThatMatrix, ImplementsCopyConstructorCorrectly)
{
  matrix<double> M1(4, 2);

  M1(1,1) = 1;M1(1,2) = 5;
  M1(2,1) = 2;M1(2,2) = 6;
  M1(3,1) = 3;M1(3,2) = 7;
  M1(4,1) = 4;M1(4,2) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test the copy constructor
  matrix<double> M2{M1};

  // Test that matrix M1 hasn't changed
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test that matrix M2 equals M1 (Copy performed correctly)
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );

  // Test that copies are independent, the default copy constructor overridden
  M2(1, 2) = 3;
  
  // Test that matrix M1 hasn't changed
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test that matrix M2 has 3 in correct position
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      if( i==1 and j==2 )
      {
        EXPECT_THAT(M2(i,j), Eq(3) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );
      }
}

TEST_F(TestThatMatrix, ImplementsCopyAssignmentCorrectly)
{
  matrix<double> M1(4, 2);

  M1(1,1) = 1;M1(1,2) = 5;
  M1(2,1) = 2;M1(2,2) = 6;
  M1(3,1) = 3;M1(3,2) = 7;
  M1(4,1) = 4;M1(4,2) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test the copy constructor
  matrix<double> M2(4,2);

  // If quick return doesn't work, them no code coverage! No Idea how else to test this.
  M1 = M1;

  M2 = M1;

  // Test that matrix M1 hasn't changed
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test that matrix M2 equals M1 (Copy performed correctly)
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );

  // Test that copies are independent, the default copy constructor overridden
  M2(2, 1) = 4;
  
  // Test that matrix M1 hasn't changed
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test that matrix M2 has 3 in correct position
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      if( i==2 and j==1 )
      {
        EXPECT_THAT(M2(i,j), Eq(4) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );
      }
}

TEST_F(TestThatMatrix, ImplementsMoveConstructorCorrectly)
{
  matrix<double> M1(4, 2);

  M1(1,1) = 1;M1(1,2) = 5;
  M1(2,1) = 2;M1(2,2) = 6;
  M1(3,1) = 3;M1(3,2) = 7;
  M1(4,1) = 4;M1(4,2) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test the Move constructor. Move is from smart pointer?
  matrix<double> M2 = move(M1);

  // Test that matrix M1 hasn't changed
  EXPECT_THAT(M1.get_cols(), Eq(0));
  EXPECT_THAT(M1.get_rows(), Eq(0));
  EXPECT_EQ(&M1(1,1), nullptr);

  // Test that matrix M2 equals old M1 (Move performed correctly)
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );

  // Test that copies are independent, the default copy constructor overridden
  M2(1, 1) = 5;
  
  // Test that matrix M2 has 5 in correct position
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      if( i==1 and j==1 )
      {
        EXPECT_THAT(M2(i,j), Eq(5) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );
      }
}

TEST_F(TestThatMatrix, ImplementsMoveAssignmentCorrectly)
{
  matrix<double> M1(4, 2);

  M1(1,1) = 1;M1(1,2) = 5;
  M1(2,1) = 2;M1(2,2) = 6;
  M1(3,1) = 3;M1(3,2) = 7;
  M1(4,1) = 4;M1(4,2) = 8;

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  // Test the Move constructor. Move is from smart pointer?
  matrix<double> M2(4, 2);

  // If quick return doesn't work, them matrix erases itself and this test cannot pass.
  M1 = move(M1);

  // Test that matrix inputed and retrieved correctly 
  for(size_t j{1}; j <= M1.get_cols(); j++)
    for(size_t i{1}; i <= M1.get_rows(); i++)
      EXPECT_THAT(M1(i,j), Eq(M1.get_rows()*(j-1) + i) );

  M2 = move(M1);

  // Test that matrix M1 hasn't changed
  EXPECT_THAT(M1.get_cols(), Eq(0));
  EXPECT_THAT(M1.get_rows(), Eq(0));
  EXPECT_EQ(&M1(1,1), nullptr);

  // Test that matrix M2 equals old M1 (Move performed correctly)
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );

  // Test that copies are independent, the default copy constructor overridden
  M2(1, 1) = 5;
  
  // Test that matrix M2 has 5 in correct position
  for(size_t j{1}; j <= M2.get_cols(); j++)
    for(size_t i{1}; i <= M2.get_rows(); i++)
      if( i==1 and j==1 )
      {
        EXPECT_THAT(M2(i,j), Eq(5) );
      }
      else
      {
        EXPECT_THAT(M2(i,j), Eq(M2.get_rows()*(j-1) + i) );
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

TEST_F(TestThatMatrix, HasDiagonalMatrixConstructorWithVector)
{
  std::vector<int> v{1, 2, 3};
  matrix<int> M(v);

  EXPECT_THAT(M.get_cols(), Eq(v.size()) );
  EXPECT_THAT(M.get_rows(), Eq(v.size()) );

  for(size_t j{1}; j <= M.get_cols(); j++)
    for(size_t i{1}; i <= M.get_rows(); i++)
      if( i==j )
      {
        EXPECT_THAT(M(i,j), Eq(v[i-1]) );
      }
      else
      {
        EXPECT_THAT(M(i,j), Eq(0) );
      }
}

TEST_F(TestThatMatrix, HasDiagonalMatrixConstructorWithInitializerList)
{
  std::vector<int> v{1, 2, 3, 4};
  matrix<int> M{1, 2, 3, 4};

  EXPECT_THAT(M.get_cols(), Eq(v.size()) );
  EXPECT_THAT(M.get_rows(), Eq(v.size()) );

  for(size_t j{1}; j <= M.get_cols(); j++)
    for(size_t i{1}; i <= M.get_rows(); i++)
      if( i==j )
      {
        EXPECT_THAT(M(i,j), Eq(v[i-1]) );
      }
      else
      {
        EXPECT_THAT(M(i,j), Eq(0) );
      }
}

TEST_F(TestThatMatrix, HasFullMatrixConstructorWithVector)
{
  matrix<int> M(2, 3, vector<int>{1, 2, 3, 4, 5, 6});

  EXPECT_THAT(M.get_cols(), Eq(3) );
  EXPECT_THAT(M.get_rows(), Eq(2) );

  EXPECT_THAT(M(1,1), Eq(1) );
  EXPECT_THAT(M(1,2), Eq(2) );
  EXPECT_THAT(M(1,3), Eq(3) );
  EXPECT_THAT(M(2,1), Eq(4) );
  EXPECT_THAT(M(2,2), Eq(5) );
  EXPECT_THAT(M(2,3), Eq(6) );
}

TEST_F(TestThatMatrix, HasFullMatrixConstructorWithInitializerList)
{
  matrix<int> M(2, 2, {1, 2, 3, 4});

  EXPECT_THAT(M.get_cols(), Eq(2) );
  EXPECT_THAT(M.get_rows(), Eq(2) );

  EXPECT_THAT(M(1,1), Eq(1) );
  EXPECT_THAT(M(1,2), Eq(2) );
  EXPECT_THAT(M(2,1), Eq(3) );
  EXPECT_THAT(M(2,2), Eq(4) );
}

// This use cases are possible
// // initializer_list constructor will be used: create a 2x2 diagonal matrix with 1, 2 on the diagonal
// cout << "matrix{1, 2}:\n";
// cout << matrix<int>{1, 2};
// // (size_t, size_t) constructor will be used: create an UNINITIALIZED 1x2 matrix
// cout << "matrix(1, 2):\n";
// cout << matrix<int>(1, 2);

