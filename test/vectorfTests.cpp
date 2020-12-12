#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cfloat>
#include <iostream>
#include <vector>

#include "vectorf.hpp"

#define DOUBLE_NEAR(x) DoubleNear((x), DBL_EPSILON)

using namespace testing;
using namespace std;

class TestThatVectorF : public Test
{
 public:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(TestThatVectorF, CanConstructAnUninitializedVectorOfSpecifiedLength)
{
  ASSERT_NO_THROW(vectorf<double> v(2));
}

TEST_F(TestThatVectorF, ConstructorSizeThrowsWheneverSizeIsNegative)
{
  ASSERT_THROW(vectorf<double> v(-1), std::bad_array_new_length);
}

TEST_F(TestThatVectorF, ConstructorSizeThrowsWheneverSizeIsZero)
{
  ASSERT_THROW(vectorf<double> v(0), std::exception);
}

TEST_F(TestThatVectorF, get_sizeWorks)
{
  vectorf<double> v(4);

  EXPECT_THAT(v.get_size(), Eq(4));
}

TEST_F(TestThatVectorF,
       CanConstructAnUninitializedVectorSepcifyingSizeANDInitialiseIt)
{
  vectorf<int> v(4);

  // Must be able to modify the element to do this
  v(1) = 1;
  v(2) = 2;
  v(3) = 3;
  v(4) = 4;

  EXPECT_THAT(v.get_size(), Eq(4));

  // Must be able to access, but not modify, the element to do this
  for (size_t i{1}; i <= v.get_size(); i++) EXPECT_THAT(v(i), Eq(i));

  // // Printing verification does work
  // std::cout << "v =\n"
  //      << v << '\n';
}

TEST_F(TestThatVectorF, ImplementsCopyConstructorCorrectly)
{
  vectorf<int> v1(4);

  v1(1) = 1;
  v1(2) = 2;
  v1(3) = 3;
  v1(4) = 4;

  // Test that vector inputed and retrieved correctly
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test the copy constructor
  vectorf<int> v2(v1);

  // Test that vector v1 hasn't changed
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test that vector v2 equals v1 (Copy performed correctly)
  for (size_t i{1}; i <= v2.get_size(); i++) EXPECT_THAT(v2(i), Eq(i));

  // Test that copies are independent, the default copy constructor overridden
  v2(1) = 3;

  // Test that vector v1 hasn't changed
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test that v2 has 3 in correct position
  for (size_t i{1}; i <= v2.get_size(); i++)
    if (i == 1)
    {
      EXPECT_THAT(v2(i), Eq(3));
    }
    else
    {
      EXPECT_THAT(v2(i), Eq(i));
    }
}

TEST_F(TestThatVectorF, ImplementsCopyAssignmentCorrectly)
{
  vectorf<int> v1(4);

  v1(1) = 1;
  v1(2) = 2;
  v1(3) = 3;
  v1(4) = 4;

  // Test that vector inputed and retrieved correctly
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test the copy constructor
  vectorf<int> v2(v1);

  // If quick return doesn't work, them no code coverage! No Idea how else to
  // test this.
  v1 = v1;

  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  v2 = v1;

  // Test that vector v1 hasn't changed
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test that vector v2 equals v1 (Copy performed correctly)
  for (size_t i{1}; i <= v2.get_size(); i++) EXPECT_THAT(v2(i), Eq(i));

  // Test that copies are independent, the default copy constructor overridden
  v2(2) = 4;

  // Test that vector v1 hasn't changed
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test that v2 has 3 in correct position
  for (size_t i{1}; i <= v2.get_size(); i++)
    if (i == 2)
    {
      EXPECT_THAT(v2(i), Eq(4));
    }
    else
    {
      EXPECT_THAT(v2(i), Eq(i));
    }
}

TEST_F(TestThatVectorF, ImplementsMoveConstructorCorrectly)
{
  vectorf<int> v1(4);

  v1(1) = 1;
  v1(2) = 2;
  v1(3) = 3;
  v1(4) = 4;

  // Test that vector inputed and retrieved correctly
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test the copy constructor
  vectorf<int> v2 = move(v1);

  // Test that vector v1 hasn't changed
  EXPECT_THAT(v1.get_size(), Eq(0));
  EXPECT_EQ(&v1(1), nullptr);

  // Test that vector v2 equals old v1 (Copy performed correctly)
  for (size_t i{1}; i <= v2.get_size(); i++) EXPECT_THAT(v2(i), Eq(i));

  // Test that copies are independent, the default copy constructor overridden
  v2(1) = 5;

  // Test that v2 has 5 in correct position
  for (size_t i{1}; i <= v2.get_size(); i++)
    if (i == 1)
    {
      EXPECT_THAT(v2(i), Eq(5));
    }
    else
    {
      EXPECT_THAT(v2(i), Eq(i));
    }
}

TEST_F(TestThatVectorF, ImplementsMoveAssignmentCorrectly)
{
  vectorf<int> v1(4);

  v1(1) = 1;
  v1(2) = 2;
  v1(3) = 3;
  v1(4) = 4;

  // Test that vector inputed and retrieved correctly
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  // Test the move assignment
  vectorf<int> v2(4);

  // If quick return doesn't work, them matrix erases itself and this test
  // cannot pass.
  v1 = move(v1);

  // Test that matrix inputed and retrieved correctly
  for (size_t i{1}; i <= v1.get_size(); i++) EXPECT_THAT(v1(i), Eq(i));

  v2 = move(v1);

  // Test that vector v1 hasn't changed
  EXPECT_THAT(v1.get_size(), Eq(0));
  EXPECT_EQ(&v1(1), nullptr);

  // Test that vector v2 equals old v1 (Copy performed correctly)
  for (size_t i{1}; i <= v2.get_size(); i++) EXPECT_THAT(v2(i), Eq(i));

  // Test that copies are independent, the default copy constructor overridden
  v2(1) = 5;

  // Test that v2 has 5 in correct position
  for (size_t i{1}; i <= v2.get_size(); i++)
    if (i == 1)
    {
      EXPECT_THAT(v2(i), Eq(5));
    }
    else
    {
      EXPECT_THAT(v2(i), Eq(i));
    }
}

TEST_F(TestThatVectorF, CannotDoCopyAssignmentWheneverSizesDonotAgree)
{
  vectorf<double> v1(4);
  vectorf<double> v2(2);

  ASSERT_THROW(v1 = v2, std::exception);
}

TEST_F(TestThatVectorF, CannotDoMoveAssignmentWheneverSizesDonotAgree)
{
  vectorf<double> v1(4);
  vectorf<double> v2(2);

  ASSERT_THROW(v1 = move(v2), std::exception);
}

TEST_F(TestThatVectorF, HasConstructorWithVector)
{
  vector<int> vec{1, 2, 3, 4, 5, 6};
  vectorf<int> v(vec);

  EXPECT_THAT(v.get_size(), Eq(6));

  EXPECT_THAT(v(1), Eq(1));
  EXPECT_THAT(v(2), Eq(2));
  EXPECT_THAT(v(3), Eq(3));
  EXPECT_THAT(v(4), Eq(4));
  EXPECT_THAT(v(5), Eq(5));
  EXPECT_THAT(v(6), Eq(6));
}

