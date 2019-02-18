function [out] = a_trim(Vt, h,gamma)
global x u gam
x(1) = Vt;
x(5) = h;
gam = gamma/57.29578;
cg = 0.25; land=1; %0 = Clean, 1= gears and flaps
u=[0.1 -10 cg land];
%name= input(�Name of Cost function file ? : �,� s�);
x(2)= 0.1; % Alpha, initial guess
x(3)= x(2) + gam; % Theta

x(4)=0; % Pitch rate set to zero
x(6)=0;
s0=[u(1) u(2) x(2)];
% Now initialize any other states and get initial cost
disp(["Initial cost = ",num2str(feval(@cost,s0))])
[s,fval]=fminsearch(@cost,s0);
x(2)=s(3); x(3)=s(3)+gam;
u(1)=s(1); u(2)=s(2);
disp(["minimum cost = ",num2str(fval)])
disp(["minimizing vector= ",num2str(s)])
temp=[length(x),length(u),x,u];
dlmwrite('test',temp);
out = 'Done';







