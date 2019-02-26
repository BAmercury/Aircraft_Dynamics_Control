function [f] = cost(state)
% Cost Function for 3-DOF Aircraft
global x u gam
u(1)= state(1);
u(2)= state(2);
x(2)= state(3);
x(3)= x(2) + gam;
time= 0.0;
[xd]=transp(time,x,u);
f= (xd(1)^2 + 100*xd(2)^2 + 10*xd(4)^2);