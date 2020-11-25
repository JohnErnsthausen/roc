# Radius of Convergence

Developer(s) Name: John Ernsthausen<br>
Date of project start: 17 September 2020

## Introduction

This software estimates the radius of convergence [(RC)](https://en.wikipedia.org/wiki/Radius_of_convergence#:~:text=The%20radius%20of%20convergence%20of%20a%20power%20series%20%C6%92%20centered,called%20the%20disk%20of%20convergence.) of a real valued power series.

## Getting started

You can build this software on a computer system equipped with the gnu c/cpp compiler. For a good
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

## Documentation for coders

The unit tests are living documentation for the functionality of this software. Read the __test-name.test-case__
as a sentence for a description of the codes functionality, and deep-dive into the test case code to understand
how to implement the functionality described.

## Documentation for users

Documentation in this repository will undergo major revisions in the coming weeks. Please bare with me
as I build the documentation for you.

## Software usage

At the moment, this software can be used as in the example __roc/examples/e_calculus.c__. The user specifies

1. The start index, kstart
2. A sufficient number of terms in the power series in TC (length at least kstart+10)
3. The scale used to produce the TCs 

The software returns the RC and the best fit line slope/intercept coefficients.

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

The code solves the linear least squares problem to best fit the power series coefficients.
The QRFactorization was implemented and used to solve the Linear Least Squares problem,
and great care was taken in the construction of the Householder transformations,
the two-norm, and the backward substitution used in solving the linear system.

This software was build somewhat for speed. However there is room for improvement. The
author anticipates calling the __roc__ subroutine thousands of times in the context
of larger software needs such as a stepsize controller in a Ordinary Differential Equation
and a Differential-Algebraic Equation solver.

This software described in [Chang and Corliss](https://dl.acm.org/doi/pdf/10.1145/355993.355995).
Linear algebra described in [Stewart](https://books.google.com/books?hl=en&lr=&id=XXzNCgAAQBAJ&oi=fnd&pg=PP1&dq=stewart+introduction+to+matrix+computations&ots=sj2tcycEub&sig=OzbhV-wuRoBZanLVSiQwtN_Gr34#v=onepage&q=stewart%20introduction%20to%20matrix%20computations&f=false) as well as __Golub, G. H., C. F. Van Loan, and Matrix Computations. "Johns Hopkins Univ." Press, Baltimore (1989)__.
Software aspects were inspired by the [MANPAK](https://www.sciencedirect.com/science/article/pii/S0898122196002040) project which I worked on with Rheinboldt.

## This repository

The folders and files for this project are as follows:

1. docs - Documentation for the project
2. refs - Reference material used for the project, including papers
3. src - Source code
4. test - Test cases
5. examples - Example use cases
