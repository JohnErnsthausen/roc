CPPSRCS = e_sevenbody.cpp \
   ../src/io.cpp \
	 ../src/qrfactorization.cpp \
	 ../src/threeterm.cpp \
	 ../src/sixterm.cpp \
	 ../src/topline.cpp \
	 ../src/roc.cpp

CPPOBJS = $(patsubst %.cpp, %.cpp.o, $(CPPSRCS))
CPPDEPS = $(patsubst %.cpp, %.cpp.d, $(CPPSRCS))

CSRCS  = ../src/dist.c \
   ../src/mathext.c \
   ../src/qrfactorization.c

COBJS  = $(patsubst %.c, %.c.o, $(CSRCS))
CDEPS  = $(patsubst %.c, %.c.d, $(CSRCS))

OBJECTS =  $(CPPOBJS) $(COBJS)
DEPENDS =  $(CPPDEPS) $(CDEPS)

WARNING := -Wall -Wextra  -Wpedantic 

CXX = g++
CXXFLAGS = -O2 -pg -std=c++17

CC = cc
CFLAGS= -O2 -pg -std=c11

LDFLAGS=-I../include
LDLIBS =-lm -lstdc++

.PHONY: all clean

all: e_sevenbody

clean:
	$(RM) $(OBJECTS) $(DEPENDS) e_sevenbody 

# Rule for the executable, depending on all the object files
e_sevenbody: $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

# We have a dependency file for each source file. Here we include them into the Makefile.
# The - here means that make will not complain if those files do not exist -- because they
# do not exist at this point.
-include $(DEPENDS)

# Generic rule to compile CPP individual source files into corresponding object file
# Demanding this rule to depend on the Makefile as well ensures that e.g. changes in
# $(CXXFLAGS) trigger a recompile as well. $< resolves to the first dependency, the
# first CPP source file. $@ is the filename of the target, here e_sevenbody
#
# We have the compiler generate those dependency files for us. For that, we extend the build
# rule. The -MMD flag generates the dependency file (%.cpp.d), which will hold in Makefile
# syntax rules making generated file (%.cpp.o) depend on the source file and any non-system
# headers it includes. That means the object file gets created automatically whenever relevant
# sources are touched. If you want to also depend on system headers, checking them for updates
# on each compile, use -MD instead.
#
# The -MP option adds empty dummy rules, which avoid errors sould header files be removed from
# the system.
#
# At the first compile run, there is no dependency information, but since the object files don't
# exist either, the compiler must run anyways. For each subsequent run, make will include the 
# auto-generated dependency files, and "do the right thing".
%.cpp.o: %.cpp 
	$(CXX) $(WARNING) $(CXXFLAGS) $(LDFLAGS) -MMD -MP -c $< -o $@

# Generic rule to compile CC individual source files into corresponding object file
# Demanding this rule to depend on the Makefile as well ensures that e.g. changes in
# $(CFLAGS) trigger a recompile as well. $< resolves to the first dependency, the
# first CC source file. $@ is the filename of the target, here e_sevenbody
%.c.o: %.c
	$(CC) $(WARNING) $(CFLAGS) $(LDFLAGS) -MMD -MP -c $< -o $@

