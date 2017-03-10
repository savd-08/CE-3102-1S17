#include <iostream>
#include <chrono>
#include "anpi.h"
#include "anpi_opt.h"
#include "anpi_ref.h"

#define ITER 300000

//using namespace anpi;

double test_opt(double x, double* terms, unsigned int term_amt){
	double res_opt = 0;
	auto begin = std::chrono::high_resolution_clock::now();
	for(int i = 0;  i < ITER; i++){
			res_opt += anpi::opt::poly_evaluator(x, terms, term_amt);

	}

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
	std::cout << "Result opt: "<< res_opt/ITER << std::endl ;
	std::cout << duration << "ns total, average : " << duration / ITER << "ns." << std::endl << std::endl;
}

double test_ref(double x, double* terms, unsigned int term_amt){
	double res = 0;
	auto begin = std::chrono::high_resolution_clock::now();
	for(int i = 0;  i < ITER; i++){
			res += anpi::ref::poly_evaluator(x, terms, term_amt);
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
	std::cout << "Result ref: "<< res/ITER << std::endl;
	std::cout << duration << "ns total, average : " << duration / ITER << "ns." << std::endl << std::endl;
}

template <typename T>
void test_ln(T center, unsigned int terms, T x){
	//logaritmo referencia
	anpi::ln_a<double> rln_a(center,terms, &anpi::ref::poly_evaluator);
	//calculo del logaritmo
	auto rbegin = std::chrono::high_resolution_clock::now();
	T rres_opt = 0;
	for(int i = 0;  i < ITER; i++){
			rres_opt += rln_a(x);
	}
	auto rend = std::chrono::high_resolution_clock::now();
	auto rduration = std::chrono::duration_cast<std::chrono::nanoseconds>(rend-rbegin).count();
	std::cout << "Reference logarithm result: "<< rres_opt/ITER << std::endl ;
	std::cout << rduration << "ns total, average : " << rduration / ITER << "ns." << std::endl << std::endl;


	//logaritmo optimizado
	anpi::ln_a<double> oln_a(center,terms, &anpi::ref::poly_evaluator);
	//calculo del logaritmo
	auto obegin = std::chrono::high_resolution_clock::now();
	T ores_opt = 0;
	for(int i = 0;  i < ITER; i++){
			ores_opt += oln_a(x);
	}
	auto oend = std::chrono::high_resolution_clock::now();
	auto oduration = std::chrono::duration_cast<std::chrono::nanoseconds>(oend-obegin).count();
	std::cout << "Optimized logarithm result: "<< ores_opt/ITER << std::endl ;
	std::cout << oduration << "ns total, average : " << oduration / ITER << "ns." << std::endl << std::endl;

}//fin de test_ln


template <typename T>
void test_cos(T center, T x, unsigned int terms){

    //Funcion coseno con Horner (referencia)
    anpi::cos_a<double> refCos(center,terms,&anpi::ref::poly_evaluator);
    T tmpresr = 0;
    //calculo del coseno
    auto rbegin = std::chrono::high_resolution_clock::now();
    for(int i = 0;  i < ITER; i++){
            tmpresr += refCos(x);
    }
    auto rend = std::chrono::high_resolution_clock::now();
    auto rduration = std::chrono::duration_cast<std::chrono::nanoseconds>(rend-rbegin).count();
    std::cout << "Reference cos result: "<< tmpresr/ITER << std::endl ;
    std::cout << rduration << "ns total, average : " << rduration / ITER << "ns." << std::endl << std::endl;

    /***************************************************************/

    //funcion coseno con Estrin (optimizada)
    anpi::cos_a<double> optCos(center,terms,&anpi::opt::poly_evaluator);
    T tmpreso = 0;
    //calculo del coseno
    auto obegin = std::chrono::high_resolution_clock::now();
    for(int i = 0;  i < ITER; i++){
            tmpreso += optCos(x);
    }
    auto oend = std::chrono::high_resolution_clock::now();
    auto oduration = std::chrono::duration_cast<std::chrono::nanoseconds>(oend-obegin).count();
    std::cout << "Optimized cos result: "<< tmpreso/ITER << std::endl ;
    std::cout << oduration << "ns total, average : " << oduration / ITER << "ns." << std::endl << std::endl;

}



int main(){
	unsigned int term_amt = 65;
	double terms[term_amt];
	double x = 2.0;
	for(int i = term_amt; i >= 0; i--){
		terms[i] = double(i);
	}
	test_ref(x, terms, term_amt);
	test_opt(x, terms, term_amt);

	double tcenter = 8;
	double tx = 10;

	test_ln(tcenter, term_amt, tx);

	/*****************************************/


	double cos_center = 0.5;
    double cos_x = 3.141592;

    test_cos(cos_center, cos_x, term_amt);



	return 0;
}
