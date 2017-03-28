%   imag vs real

load laguerre.txt
plot (laguerre(:,1),laguerre(:,2),";Laguerre's method;","color",'r')
hold on
load muller.txt
plot (muller(:,1),muller(:,2),";Muller's method;","color",'g')
title ('ROOTS')
xlabel ('Real')
ylabel ('Imaginary')

pause(45)

