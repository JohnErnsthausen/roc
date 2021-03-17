#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "exceptions.hpp"

using namespace testing;

TEST(TestExceptions, HandleWhatMethodOfBaseClassException)
{
  try
  {
    std::string message{"I am a message"};
    throw sayMessage(message);
  }
  catch (const sayMessage& e)
  {
    ASSERT_STREQ(e.what(), "I am a message");
  }
}
