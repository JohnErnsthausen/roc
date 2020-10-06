cmake_minimum_required(VERSION 3.10)

# Radius Of Convergence
project(ROC)

# You can tweak some common (for all subprojects) stuff here. For example:

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif ()

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_COLOR_MAKEFILE   ON)

# Remove 'lib' prefix for shared libraries on Windows
if (WIN32)
  set(CMAKE_SHARED_LIBRARY_PREFIX "")
endif ()

set(CMAKE_BUILD_PARALLEL_LEVEL 8)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")

set(PROJECT_LINK_LIBS libroc.a)
message("-- Using this project source dir: [${PROJECT_SOURCE_DIR}]")
link_directories(${PROJECT_SOURCE_DIR}/build)

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB ROC_SRC src/*.c src/*.cpp)
add_library(roc STATIC ${ROC_SRC})

find_package(GTest REQUIRED)
if(GTEST_FOUND)
  set(Gtest_FOUND TRUE)
endif()
if(GTest_FOUND)
  include(GoogleTest)
  enable_testing()
  add_executable(run test/MainTests.cpp test/HouseholderTests.cpp )
  target_link_libraries(run GTest::GTest)
  target_link_libraries(run ${PROJECT_LINK_LIBS})
  gtest_discover_tests(run)
endif()
