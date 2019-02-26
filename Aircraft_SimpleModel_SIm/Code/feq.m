%%feq
function [F] = feq(x)
global num_1 num_2 den_1 den_2 den_3
global cdcls cl0 cla


dd = cdcls*cl0*cla;

F = tan(x*pi/180)*(den_1+den_2*x+den_3*x^2)+(num_1+num_2*x);
F = F^2;
end