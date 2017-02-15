#include "tarea_1.h"

int main(int argc, char** argv){

	cout << "\n--- Tarea 1 - Análisis Numérico ---\n\n";
	cout << "------------Ejercicio 2------------\n";
	cout << "Valor Real: 1\n";
	cout << "Precisión simple: " << ej_2(100000, 0.00001f) << "\n";
	cout << "Precisión doble:  " << ej_2(100000, 0.00001)  << "\n\n";

	quad_ec_sol<float> res_float = ej_3(1.0f, 80.0f, 16.0f);
	quad_ec_sol<double> res_double = ej_3(1.0, 80.0, 16.0);

	cout << "------------Ejercicio 3------------\n";
	cout << "Parámetros de prueba: a = 1; b = 80; c = 16; (80²) >> 4*1*16 \n";
	cout << "Valor Real: x1 = -0,02005025157; x2 = -79.79949748\n";
	cout << "Precisión simple: x1 = " << res_float.x1 << "; x2 = " << res_float.x2 << "\n";
	cout << "Precisión doble:  x1 = " << res_double.x1 << "; x2 = " << res_double.x2 << "\n\n";

	return 0;
}

template <typename fp>
fp ej_2 (int n, fp x) {

	//Se inicializa un acumulador en cero
	fp acum = fp(0);		

	//Se incrementa el acumulador en x, n veces 
	for(int i = 0; i < n; i++){		
		acum += x;			
	}

	return acum;
}


template <typename fp>
quad_ec_sol<fp> ej_3 (fp a, fp b, fp c) {

	//Estructura que contiene las soluciones (x1, x2)
	quad_ec_sol<fp> result; 

	//Se obtiene el discriminante (delta) de la ecuación
	//sqrt(b² - 4*a*c)
	fp delta = pow(b, fp(2)); 
	delta -= fp(4) * a * c;
	delta = sqrt(delta);

	//Se calcula de antemano el denominador de las soluciones
	//2*a
	fp denom = fp(2)*a; 

	//Se calcula una de las soluciones
	//(-b + delta)/(2*a)
	fp x1 = -b + delta;
	x1 = x1/denom;
	
	//Se calcula la solución restante
	//(-b - delta)/(2*a)
	fp x2 = -b - delta;
	x2 = x2/denom;

	//Se colocan los resultados en la estructura para tal fin
	result.x1 = x1;
	result.x2 = x2;	

	return result;
}