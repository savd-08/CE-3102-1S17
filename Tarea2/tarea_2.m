#Arreglo de exponente para la h
exponente = 0:20;
#Arreglo para el calculo de las h
h_arreglo = 10 .^ (-exponente);
#Valor de donde se evalua la funcion
x = 1;
#Valor real de la derivada
valor_real = 0.9;


#Funcion que calcula el valor de la funcion en cuestion en una posicion x
function retval = f_x  (x)
  if(isscalar(x))
    retval = (0.3*(x**4)) - (0.15*(x**2));
  else
    error("f_x: se espera un argumento escalar");
  endif
endfunction
#Funcion que calcula el valor por diferencia centrada
function retval = d_centrada (x,h)
  retval = (f_x(x+h) - f_x(x-h))/(2*h);
endfunction
#funcion que calcula el valor de la diferencia hacia atras
function retval = d_atras (x,h)
  retval = (f_x(x) - f_x(x-h))/h;
endfunction
#funcion que calcula el valor de la diferencia hacia delante
function retval = d_adelante (x,h)
  retval = (f_x(x+h) - f_x(x))/h;
endfunction

#calulo de los errores
err_centrada = abs(valor_real .- arrayfun(@d_centrada, x, h_arreglo));

err_atras = abs(valor_real .- arrayfun(@d_atras, x, h_arreglo));

err_adelante = abs(valor_real .- arrayfun(@d_adelante, x, h_arreglo));

#graficacion de los erroes respecto a la diferencia
loglog(h_arreglo, err_centrada);
xlabel("log de la diferencia");
ylabel("log del error");
title("Diferencias Centradas")

loglog(h_arreglo, err_atras);
xlabel("log de la diferencia");
ylabel("log del error");
title("Diferencias Hacia Atras")

loglog(h_arreglo, err_adelante);
xlabel("log de la diferencia");
ylabel("log del error");
title("Diferencias Delante")

