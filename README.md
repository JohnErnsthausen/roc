# Radius of Convergence

Developer(s) Name: John Ernsthausen<br>
Date of project start: 17 September 2020

## This repo uses

[![Build Status](https://travis-ci.org/JohnErnsthausen/roc.svg?branch=master)](https://travis-ci.org/JohnErnsthausen/roc)
[![Coverage Status](https://coveralls.io/repos/github/JohnErnsthausen/roc/badge.svg?branch=master)](https://coveralls.io/github/JohnErnsthausen/roc?branch=master)

[Coverage from plain gcov](http://www.johnernsthausen.com/experiences/coverage/)

## Introduction

This software estimates the radius of convergence [(RC)](https://en.wikipedia.org/wiki/Radius_of_convergence#:~:text=The%20radius%20of%20convergence%20of%20a%20power%20series%20%C6%92%20centered,called%20the%20disk%20of%20convergence.) of a real valued power series.

## Getting started

You can build this software on a computer system equipped with the gnu cpp compiler. For a good
tutorial on equipping your system with this compiler, see
[Barak Shoshany](http://baraksh.com/CSE701/notes.php#visual-studio-code).

For convince and to run the test suite, additionally equip your system with

1. [git](https://git-scm.com/downloads)
1. [cmake](https://cmake.org/download)
1. the Google unit testing suite [GTest and GMock](https://github.com/google/googletest)

Properly integrate linting, code coverage, profiling, and performance testing by
additionally equipping your system with the following development tools

1. clang-check, clang-format, and clang-tidy aka [ClangTools](http://clang.llvm.org/docs/ClangTools.html)
1. [valgrind](https://valgrind.org)
1. gprof aka [Performance Analysis tools](https://en.wikipedia.org/wiki/List_of_performance_analysis_tools)
1. [gcov](https://en.wikipedia.org/wiki/Gcov)

We implement continuous integration and build testing through TravisCI.

A Ruby Rakefile is offered for convenience for those users who can take advantage of it. Observe its
functionality through the command `rake -T`.

An IDE such as VSCode will provide many tools for convenience. Editors such as GVim
provide an alternative to the perks and advantages offered through VSCode.

With the software described above in place, there are two ways to compile the code. The standard way assumes
the developer is in the top level directory roc/. Then

> mkdir build

> cd build

> cmake ..

> make

The Ruby shortcut to do all this in one simple command is

> rake make

The Ruby Rakefile has a lot of cool functionality. In the top level directory, simply type at the command
line

> rake -T

to learn all the advantage of this small but mighty script. For example, format the code with

> rake clang

or run all the unit tests with

> rake test

> rake run

Clean up the repository with

> rake clean

> rake clobber

Examples reside in the directory roc/examples. Change to this directory and type

> rake compile[e_calculus]

and a Makefile will be constructed and executed. To execute the software built for this example type

> ./e_calculus

Clean up the directories with a

> make clean

> rake clean

To remove the Makefile, type

> rake clobber

## Example Calculus

We present an example form [a calculus tutorial](https://tutorial.math.lamar.edu/Classes/CalcII/PowerSeries.aspx)
and analyse

```
f(x) = sum_0^{infty} \frac{2^n}{n} (4 x -8)^n on 15/8 <= x < 17/8 has Rc=1/8.
```

From roc directory

1. cd examples
2. rake compile[e_calculus]
3. ./e_calculus
4. rake clean

## Example DAE

We present an example with Taylor series data from the solution of a Differential Algebraic Equation. From roc directory 

1. cd examples
2. rake compile[e_laynewatson]
3. ./e_laynewatson
4. gprof e_laynewatson gmon.out > gprof_output.txt
4. rake clean

## Documentation for coders

The unit tests are living documentation for the functionality of this software. Read the __test-name.test-case__
as a sentence for a description of the codes functionality, and deep-dive into the test case code to understand
how to implement the functionality described.

## Documentation for users

Documentation in this repository is available in the `docs` directory.

## Obtain a Taylor Series of a function 

Download [FADBAD++ 2.1](http://www.fadbad.com/fadbad.html)

In the program GetTS.cpp in the top level directory alter the following to your function

```bash
T<double> func(const T<double>& x)
{
  T<double> y = 1.0 - x;
  return 1.0/y;
}
```

```bash
double x0{0.3}, hbar{1.0};
```

Put this program into FADBAD++/examples and compile with

> gcc -I.. GetTS.cpp -o ex.exe -lstdc++ -lm

and run to obtain a 30 term Taylor Series output to STDOUT

> ./ex.exe


## Algorithms used

The code uses LAPACK/DGELS and LAPACK/dnrm2.

This software was built for speed. The
author anticipates calling the __roc__ subroutine thousands of times in the context
of larger software needs such as a stepsize controller in a Ordinary Differential Equation
and a Differential-Algebraic Equation solver.

This software described in [Chang and Corliss](https://dl.acm.org/doi/pdf/10.1145/355993.355995).

## This repository

The folders and files for this project are as follows:

1. docs - Documentation for the project
2. refs - Reference material used for the project, including papers
3. cmake - Cmake build resources used in CMakeLists.txt
4. src - Source code
5. include - Include files
6. test - Test cases
7. examples - Example use cases. See the Rakefile.rb for building make files to compile examples
8. GetTS.cpp - Compute Taylor coefficients of a function
9. LICENSE - License information
10. Rakefile.rb - Ruby tasks for your convenience, although not necessary to use them
11. README.md - This readme file

