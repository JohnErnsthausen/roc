cmake_minimum_required(VERSION 3.16)

# Radius Of Convergence
project(ROC)

# You can tweak some common (for all subprojects) stuff here. For example:

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
option(BUILD_TESTS "Build test executable" ON)
option(ACCEPTANCE_TESTS "Build acceptance test executable" OFF)
option(ENABLE_COVERAGE "Enable code coverage" OFF)

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

if(ENABLE_COVERAGE AND NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT BUILD_TESTS)
    message(FATAL_ERROR "Tests must be enabled for code coverage!")
  endif(NOT BUILD_TESTS)

  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake") 
  include(CodeCoverage)

  append_coverage_compiler_flags()
  set(COVERAGE_EXCLUDES
      ${CMAKE_SOURCE_DIR}/docs/\\*
      ${CMAKE_SOURCE_DIR}/examples/\\*
      ${CMAKE_SOURCE_DIR}/refs/\\*
      /usr/\\*)

  setup_target_for_coverage(NAME coverage EXECUTABLE ctest DEPENDENCIES
                            coverage)
else(ENABLE_COVERAGE AND NOT CMAKE_CONFIGURATION_TYPES)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL
                                              "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /WX")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
  endif(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL
                                                 "GNU")
endif(ENABLE_COVERAGE AND NOT CMAKE_CONFIGURATION_TYPES)

if(BUILD_TESTS)
  find_package(LAPACK)
  if(NOT LAPACK_LIBRARIES)
    message(FATAL_ERROR "Cannot find LAPACK libraries")
  endif()
  find_package(BLAS)
  if(NOT BLAS_LIBRARIES)
    message(FATAL_ERROR "Cannot find BLAS libraries")
  endif()
  find_package(GTest REQUIRED)
  if(GTEST_FOUND)
    set(GTest_FOUND TRUE)
  endif()
  if(GTest_FOUND)
    enable_testing()

    file(GLOB ROC_TESTS "${PROJECT_SOURCE_DIR}/test/*.cpp")
    set(TEST_EXE_NAME run)
    add_executable(${TEST_EXE_NAME} ${ROC_TESTS})
    target_include_directories(${TEST_EXE_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /home/john/ipopt/include)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /usr/lib/gcc/x86_64-linux-gnu/9)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /usr/lib/gcc/x86_64-linux-gnu/9/../../../x86_64-linux-gnu)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /usr/lib/gcc/x86_64-linux-gnu/9/../../../../lib)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /lib/x86_64-linux-gnu)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /lib/../lib)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /usr/lib/x86_64-linux-gnu)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /usr/lib/../lib)
    target_include_directories(${TEST_EXE_NAME} PUBLIC /usr/lib/gcc/x86_64-linux-gnu/9/../../..)
    target_link_libraries(${TEST_EXE_NAME} GTest::GTest)
    target_link_libraries(${TEST_EXE_NAME} roc)
    target_link_libraries(${TEST_EXE_NAME} /home/john/ipopt/lib/libipopt.a)
    target_link_libraries(${TEST_EXE_NAME} lapack)
    target_link_libraries(${TEST_EXE_NAME} blas)
    target_link_libraries(${TEST_EXE_NAME} gfortran)
    target_link_libraries(${TEST_EXE_NAME} m)
    target_link_libraries(${TEST_EXE_NAME} gcc_s)
    target_link_libraries(${TEST_EXE_NAME} quadmath)
    gtest_discover_tests(${TEST_EXE_NAME})
  endif()
endif(BUILD_TESTS)

if(ACCEPTANCE_TESTS)
  find_package(LAPACK)
  if(NOT LAPACK_LIBRARIES)
    message(FATAL_ERROR "Cannot find LAPACK libraries")
  endif()
  find_package(BLAS)
  if(NOT BLAS_LIBRARIES)
    message(FATAL_ERROR "Cannot find BLAS libraries")
  endif()
  find_package(GTest REQUIRED)
  if(GTEST_FOUND)
    set(GTest_FOUND TRUE)
  endif()
  if(GTest_FOUND)
    enable_testing()

    file(GLOB ACCEPTANCE_TESTS "${PROJECT_SOURCE_DIR}/accept/*.cpp")
    set(TEST_EXE_NAME accept)
    add_executable(${TEST_EXE_NAME} ${ACCEPTANCE_TESTS})
    target_include_directories(${TEST_EXE_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${TEST_EXE_NAME} GTest::GTest)
    target_link_libraries(${TEST_EXE_NAME} roc)
    target_link_libraries(${TEST_EXE_NAME} lapack)
    target_link_libraries(${TEST_EXE_NAME} blas)
    gtest_discover_tests(${TEST_EXE_NAME})
  endif()
endif(ACCEPTANCE_TESTS)

