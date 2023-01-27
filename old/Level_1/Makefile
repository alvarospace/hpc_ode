include /gpfs/projects/bsc21/bsc21808/Cantera-Alya/cantera/INTEL/include/cantera/Cantera.mak

CC=icc
CXX=icpc
RM=rm -f
CXXFLAGS=-vec-report0 -diag-disable 1478 -pthread -O3 -g -std=c++0x -openmp
CPPFLAGS=$(CANTERA_INCLUDES)
LDFLAGS=
LDLIBS=$(CANTERA_LIBS)

SRCS=serial.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: serial

clean:
	$(RM) $(OBJS) serial

dist-clean: clean
	$(RM) *~ 
