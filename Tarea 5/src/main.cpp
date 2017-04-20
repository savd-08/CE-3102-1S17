#include "qr_decomp.hpp"

int main(int argc, char *argv[]){
	std::cout.precision(std::numeric_limits<double>::digits10 + 1);

	double a1[] = {1,2,2,1};								//Sol: {1/3,1/3} 
	double a2[] = {3.0,2.0,1.0,5.0,3.0,4.0,1.0,1.0,-1.0};	//Sol: {-4,6,1} 
	double a3[] = {6,3,-1,1,3,1,1,1,1,2,1,0,2,3,4,0};		//Sol: {3,-4,1,2} 
	double a4[] = {1,2,3,4,5,6,7,8,9}; 						//Singular
	double a5[] = {0.96,68,0.69,69,0.24,56,0.18,52,0.58,95,0.71,92,0.90,107,0.81,142};//Mal condicionada

	anpi::Matrix<double> A(2, 2, a1); //Variar el vector de datos y las dimensiones
	std::vector<double> x(2); //Cambiar el tamaño de acuerdo con dimensiones de la matriz

	//Vectores b para distintas matrices A
	std::vector<double> b1 = {1, 1};
	std::vector<double> b2 = {1.0,2.0,1.0};
	std::vector<double> b3 = {3,4,-4,-2};
	std::vector<double> b4 = {10,11,12};
	std::vector<double> b5 = {4,0.11,25,0.5};

	//------------------------------------------PRUEBAS QR------------------------------------------
	std::cout << "Descomposición QR" << std::endl;
    anpi::Matrix<double> Q;
	anpi::Matrix<double> R;

	//******Prueba de descomposición******
	anpi::qr(A, Q, R);	

	anpi::print_Mat(A, 'A');
	std::cout << std::endl;
	anpi::print_Mat(Q, 'Q');
	std::cout << std::endl;
	anpi::print_Mat(R, 'R');
	std::cout << std::endl;

	//Comparación de Ax y Q*R
	double norm = anpi::testQR(A, Q, R);
	std::cout << "||A-Q*R||² = " << norm << std::endl;

	//*****Resolución de sistemas de ecuaciones*****
	std::cout << std::endl << "A*x = b" << std::endl << std::endl;
	anpi::print_Mat(A, 'A');
	std::cout << std::endl;
	anpi::print_vec(b1, 'b');//Cambiar ID de b
	std::cout << std::endl; 
	
	anpi::solveQR(A, x, b1);//Cambiar ID de b
	anpi::print_vec(x, 'x');


    return 0;
}