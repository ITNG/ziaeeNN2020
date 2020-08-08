#!/bash/bin

CXX = g++
CXXFLAGS = -std=c++11 -O2 -lm -lgsl -lgslcblas -fopenmp # -lfftw3 -Wall -Wextra -g
 

prog: main.o lib.o  sdde_solver.o
	$(CXX)  -o $@  $^ $(CXXFLAGS)

main.o : main.cpp lib.h sdde_solver.h
	$(CXX) -c main.cpp  $(CXXFLAGS)

lib.o: lib.cpp lib.h sdde_solver.h 
	$(CXX) -c lib.cpp $(CXXFLAGS)

sdde_solver.o: sdde_solver.cpp sdde_solver.h lib.h
	$(CXX) -c sdde_solver.cpp $(CXXFLAGS)	
	
.PHONEY: clean
clean:
	rm -rf  *.pyc ../data/text/* ../data/fig/*
	rm -f  *.out *.pyc ../data/npz/*

.PHONEY: eradicate
eradicate: clean
	rm -rf *.o  prog
	# rm -f networks/*.txt

.PHONY: clean, eradicate
