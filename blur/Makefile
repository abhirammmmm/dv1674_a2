# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -O3 -Wall -Wunused
PTHREADS=-pthread

all: blur blur_par

blur: matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur

filters: matrix filters.hpp filters.cpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o

matrix: matrix.hpp matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ppm: ppm.hpp ppm.cpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

filters_par.o: filters_par.cpp filters_par.hpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c filters_par.cpp -o filters_par.o

blur_par: blur_par.cpp matrix.o ppm.o filters_par.o
	$(CXX) $(CXXFLAGS) blur_par.cpp matrix.o ppm.o filters_par.o -o blur_par $(PTHREADS)

clean:
	rm -rf blur blur_par *.ppm *.o *.dSYM 2> /dev/null
