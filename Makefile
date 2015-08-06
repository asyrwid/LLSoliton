all:
	c++ main.cpp permutator.cpp --std=c++11 -o permutator.out -O3 --fast-math -fopenmp
run: all
	./permutator.out < first_number 1>Output 

