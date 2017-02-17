function [trueval, result, Et, Ea, terms] = euler( n , x )
  # se establece una precision mayor que la de los numeros significativos
  presicion = output_precision();
  # se obtiene el valor ¨real¨, tomado de la constante del lenguaje
  trueval = e ^ x;
  # factorial
  fact = 1;
  # resultado
  result = 1;
  # resultado anterior
  preresult = -1;
  # valor que se espera que tenga con los numero significativos
  valexpect = floor(trueval * (10^n));
  # numero de terminos que se usan para llegar al valor anterior
  terms = 1;
  # temporal para determinar si se llego al valor especificado
  temp = 0;
  
  if(n > presicion)
    disp("la cantidad de numeros significativos son mayores a la precision");
    trueval = -1;
    result = -1;
    Et=-1;
    Ea=-1;
    terms=-1;
    return;
  endif
  
  do
 
    preresult = result;
    
    # agrega un termino al resultado
    result = result + ( (x^terms) / fact );
    
    # itera los valores par el proximo termino
    terms = terms + 1;
    
    fact = fact * terms;
    
    temp = floor(result * (10^n));
    
    # calcula error aproximado de la iteracion
    Ea = (1 - (preresult/result))*100
  
  until temp == valexpect
  # error verdadero
  Et = (1 - (result/trueval))*100
  
endfunction

output_precision(10);  