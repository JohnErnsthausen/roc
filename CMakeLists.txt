cmake_minimum_required(VERSION 3.16)

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
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")

message("-- Using this project source dir: [${PROJECT_SOURCE_DIR}]")

file(GLOB all_SRCS
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.c"
        )
add_library(roc SHARED ${all_SRCS})
target_include_directories(roc PRIVATE ${PROJECT_SOURCE_DIR}/include)

find_package(GTest REQUIRED)
if(GTEST_FOUND)
  set(Gtest_FOUND TRUE)
endif()
if(GTest_FOUND)
  file(GLOB ROC_TESTS "${PROJECT_SOURCE_DIR}/test/*.cpp")
  set(TEST_EXE_NAME run)
  enable_testing()
  include(GoogleTest)
  add_executable(${TEST_EXE_NAME} ${ROC_TESTS})
  target_include_directories(${TEST_EXE_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(${TEST_EXE_NAME} GTest::GTest)
  target_link_libraries(${TEST_EXE_NAME} roc)
  gtest_discover_tests(${TEST_EXE_NAME})
endif()
