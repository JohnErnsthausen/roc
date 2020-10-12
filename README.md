# Radius of Convergence

Developer(s) Name: John Ernsthausen<br>
Date of project start: 17 September 2020

## Introduction

Given a power series, this software computes its radius of convergence
[(RC)](https://en.wikipedia.org/wiki/Radius_of_convergence#:~:text=The%20radius%20of%20convergence%20of%20a%20power%20series%20%C6%92%20centered,called%20the%20disk%20of%20convergence.). This software is written in c/cpp.

## Getting started

Get started by installing the gnu c/cpp compiler on your system. While the following development tools
are not required, they make development and evaluation easier. Software development is facilitated with
the Ruby scripting language, make, cmake, clang-format, valgrind, git, and the Google unit testing suite
GTest and GMock. An IDE such as VSCode will provide many tools for convenience. For the hard core developer,
editors such as GVim provide an alternative to the perks and advantages of an IDE.

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

## Algorithms used

The software was built for numerical stability. Great care was taken to solve the linear
least squares problem to best fit the power series coefficients. The QRFactorization
was implemented and used to solve the Linear Least Squares problem, and great care for
numerical stability was taken in the construction of the Householder transformations,
the two-norm, and the backward substitution used in solving the linear system.

This software was build somewhat for speed. However there is room for improvement. The
author anticipates calling the __roc__ subroutine thousands of times in the context
of larger software needs such as a stepsize controller in a Ordinary Differential Equation
and a Differential-Algebraic Equation solver.

This software described in [Chang and Corliss](https://dl.acm.org/doi/pdf/10.1145/355993.355995).
Linear algebra described in [Stewart](https://books.google.com/books?hl=en&lr=&id=XXzNCgAAQBAJ&oi=fnd&pg=PP1&dq=stewart+introduction+to+matrix+computations&ots=sj2tcycEub&sig=OzbhV-wuRoBZanLVSiQwtN_Gr34#v=onepage&q=stewart%20introduction%20to%20matrix%20computations&f=false) as well as __Golub, G. H., C. F. Van Loan, and Matrix Computations. "Johns Hopkins Univ." Press, Baltimore (1989)__.
Software aspects were inspired by the [MANPAK](https://www.sciencedirect.com/science/article/pii/S0898122196002040) project I worked on Rheinboldt.

## This repository

The folders and files for this project are as follows:

docs - Documentation for the project
refs - Reference material used for the project, including papers
src - Source code
test - Test cases
examples - Example use cases
