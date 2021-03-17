CPPSRCS =
CPPOBJS = $(patsubst %.cpp, %.cpp.o, $(CPPSRCS))
CPPDEPS = $(patsubst %.cpp, %.cpp.d, $(CPPSRCS))

CSRCS  = e_calculus.c \
   ../src/dist.c \
   ../src/mathext.c \
   ../src/qrfactorization.c \
	 ../src/roc.c

COBJS  = $(patsubst %.c, %.c.o, $(CSRCS))
CDEPS  = $(patsubst %.c, %.c.d, $(CSRCS))

OBJECTS =  $(COBJS)
DEPENDS =  $(CDEPS)

WARNING := -Wall -Wextra  -Wpedantic 

CXX = g++
CXXFLAGS =$(RESEARCH) -g --std=c++11

CC = cc
CFLAGS=-g -std=c11

LDFLAGS=-I../include
LDLIBS =-lm

.PHONY: all clean

all: e_calculus

clean:
	$(RM) $(OBJECTS) $(DEPENDS) e_calculus 

e_calculus: $(OBJECTS)
	$(CC) $(WARNING) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

-include $(DEPENDS)

%.c.o: %.c
	$(CC) $(WARNING) $(CFLAGS) $(LDFLAGS) -MMD -MP -c $< -o $@
