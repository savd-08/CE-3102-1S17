#include <iostream>
#include <chrono>
#include "anpi.h"

#define ITER 250000

//using namespace anpi;

int main(){
	//double terms[10] = {11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0};
	double terms[3] = {4.0, 3.0, 2.0};
	double res = 0;
	auto begin = std::chrono::high_resolution_clock::now();
	for(int i = 0;  i < ITER; i++){
			res += anpi::opt::poly_evaluator(3.0, terms, 3);
	}
	
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
	std::cout << duration << "ns total, average : " << duration / ITER << "ns.\n" << std::endl;
	std::cout << "Result: "<< res/ITER << std::endl;
	return 0;
}