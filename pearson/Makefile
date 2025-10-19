# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -Wunused -Wall -Wextra

all: pearson verify

pearson: vector dataset analysis pearson.cpp 
	$(CXX) $(CXXFLAGS) pearson.cpp vector.o dataset.o analysis.o -o pearson

pearson_par: analysis.o dataset.o vector.o pearson_par.o
	$(CXX) $(CXXFLAGS) pearson_par.o analysis.o dataset.o vector.o -o pearson_par -lpthread

pearson_par.o: pearson_par.cpp
	$(CXX) $(CXXFLAGS) -c pearson_par.cpp -o pearson_par.o

analysis: vector analysis.hpp analysis.cpp
	$(CXX) $(CXXFLAGS) -c analysis.cpp -o analysis.o

dataset: vector dataset.hpp dataset.cpp
	$(CXX) $(CXXFLAGS) -c dataset.cpp -o dataset.o

vector: vector.hpp vector.cpp
	$(CXX) $(CXXFLAGS) -c vector.cpp -o vector.o

verify: verify.c
	$(CC) verify.c -o verify

clean:
	rm -rf verify pearson *.o *.dSYM 2> /dev/null
