#include <boost/program_options.hpp>
#include <cmath>
#include <chrono>
#include <iostream>
#include "liebmann.hpp"
#include <limits>
#include <string>

int main(int argc, char** argv) {

  //Valor constante para NaN
  static const double nan = sqrt(std::numeric_limits<double>::quiet_NaN());

  //Instanciar las variables correspondientes para los argumentos
  //No se inicializan debido a que program_options asigna valores por defecto;
  int plateSize;
  bool heatFlux;
  bool optimize;
  bool ui;
  double lambda;
  double thermCond;
  double relativeError;

  //Temperaturas tienen valor NaN para representar estado aislado
  double leftT = nan;
  double rightT = nan;
  double topT = nan;
  double bottomT = nan;

  //Cuenta de lados aislados
  uint8_t isolatedCount = 4;

  try{

    namespace po = boost::program_options;

    //Se declaran las opciones del programa
    po::options_description desc("Program options");
    desc.add_options()
        ("help,h", "print help message")
        ("optimize,o", po::bool_switch(&optimize)->default_value(false), "optimize execution time using OpenMP")
        ("heat-flux,v", po::bool_switch(&heatFlux)->default_value(false), "calculate heat flux vector field")
        ("user-interface,u", po::bool_switch(&ui)->default_value(false), "launch user interface")
        ("plate-size,s", po::value<int>(&plateSize)->default_value(10), "square plate side size (cm)")
        ("lambda,y", po::value<double>(&lambda)->default_value(1.5), "relaxation factor (1<y<2)")
        ("thermal-coeff,k", po::value<double>(&thermCond)->default_value(0.5), "thermal conductivity coefficient")
        ("relative-error,e", po::value<double>(&relativeError)->default_value(0.5), "relative error termination condition")
        ("t-left,l", po::value<double>(&leftT), "left border temperature")
        ("t-right,r", po::value<double>(&rightT), "right border temperature")
        ("t-top,t", po::value<double>(&topT), "top border temperature")
        ("t-bottom,b", po::value<double>(&bottomT), "bottom border temperature")
    ;
    po::variables_map vm;


    try {
      //Se parsean las descripciones de las opciones
      po::store(po::parse_command_line(argc, argv, desc), vm);

      //Imprime las descripciones
      if (vm.count("help")) {
          std::cout << desc << std::endl;
          return 0;
      }

      po::notify(vm);

    } catch(po::error& e){
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      std::cerr << desc << std::endl;
      return 1;
    }

    //Verificar que se hayan establecido al menos 2 temperaturas
    if(vm.count("t-left")){
      isolatedCount -= 1;
    }
    if(vm.count("t-right")){
      isolatedCount -= 1;
    }
    if(vm.count("t-top")){
      isolatedCount -= 1;
    }
    if(vm.count("t-bottom")){
      isolatedCount -= 1;
    }

    if (isolatedCount > 2)
    {
      std::cout << "There cannot be more than 2 isolated borders." << std::endl;
      return 0;
    }

    if(vm.count("lambda")){
      if(lambda < 1 || lambda > 2){
        std::cout << "Relaxation factor value must be bewtween 1 and 2." << std::endl;
        return 0;
      }
    }

  } catch(std::exception& e) {
    std::cerr << "Unhandled Exception reached the top of main: "
              << e.what() << ", application will now exit" << std::endl;
    return 1;

  }

  auto begin = std::chrono::high_resolution_clock::now();


  //Matriz para método de Liebmann
  anpi::Matrix<double> M(plateSize,plateSize,0.0);
  anpi::Matrix< std::pair<double,double> > M2(plateSize,plateSize,std::pair<double,double>(0.0,0.0));



  //Ejecución de método de Liebmann, con opciones ingresadas por el usuario
  anpi::liebmann(M, topT, rightT, bottomT, leftT, lambda, relativeError, optimize);

  /*for (size_t i = 0; i < plateSize; i++) {
    for (size_t j = 0; j < plateSize; j++) {
      std::cout << M(i,j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;*/


  //Se calculan los vectores de flujo de calor si el usuario lo desea
  if(heatFlux){
    anpi::heatflux(M, M2, topT, rightT, bottomT, leftT, thermCond, optimize);
    std::cout << "Heat Flux." << std::endl;
    /*for (size_t i = 0; i < plateSize; i++) {
      for (size_t j = 0; j < plateSize; j++) {
        std::cout << "(" << M2(i,j).first << "," << M2(i,j).second << ")" << " ";
      }
      std::cout << std::endl;
    }*/
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "Execution time: " << duration << "ns." << std::endl << std::endl;

  //escritura de archivos
  anpi::writeHeatMap(M);
  anpi::writeFlux(M2);

  //Interfaz gráfica
  if(ui){
  	int sys_msg = system("cd ui && python flux_heatmap.py");
	if(sys_msg == -1){
		std::cout << "Command line failed" << std::endl;
	}
  }

  return 0;
}
