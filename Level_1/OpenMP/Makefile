include /home/alumno.upv.es/almousa/install/cantera-2.5.1-sundials-5.3.0/include/cantera/Cantera.mak

CC=mpicc
CXX=mpicxx
RM=rm -f
CXXFLAGS=-fopenmp -g
CPPFLAGS=$(CANTERA_INCLUDES)
LDFLAGS=
LDLIBS=$(CANTERA_LIBS)

SRCS=main.cpp utils.hpp
OBJS=$(subst .cpp,.hpp,.o,$(SRCS))

all: main

clean:
	$(RM) $(OBJS) main

dist-clean: clean
	$(RM) *~ 
