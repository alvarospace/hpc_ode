include /home/almousa/install/cantera_alone/include/cantera/Cantera.mak

CC=gcc
CXX=g++
RM=rm -f
CXXFLAGS=-fopenmp -g
CPPFLAGS=$(CANTERA_INCLUDES)
LDFLAGS=
LDLIBS=$(CANTERA_LIBS) -lyaml-cpp

SRCS=main.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: main

main: main.cpp utils.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) main.cpp $(LDLIBS) -o main

clean:
	$(RM) main

dist-clean: clean
	$(RM) *~ 
