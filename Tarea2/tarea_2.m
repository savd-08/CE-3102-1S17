h_arr = [1, 1e-1 , 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10];

x = 1;

valor_real = 0.9;



function retval = f_x  (x)
  if(isscalar(x))
    retval = (0.3*(x**4)) - (0.15*(x**2));
  else
    error("f_x: se espera un argumento escalar");
  endif
endfunction

function retval = d_centrada (x,h)
  retval = (f_x(x+h) - f_x(x-h))/(2*h);
endfunction

function retval = d_atras (x,h)
  retval = (f_x(x) - f_x(x-h))/h;
endfunction

function retval = d_adelante (x,h)
  retval = (f_x(x+h) - f_x(x))/h;
endfunction


err_centrada = valor_real .- arrayfun(@d_centrada, x, h_arr);

err_atras = valor_real .- arrayfun(@d_atras, x, h_arr);

err_adelante = valor_real .- arrayfun(@d_adelante, x, h_arr);


disp("Error centrada: "), disp(err_centrada);

disp("Error atras: "), disp(err_atras);

disp("Error adelante: "), disp(err_adelante);

