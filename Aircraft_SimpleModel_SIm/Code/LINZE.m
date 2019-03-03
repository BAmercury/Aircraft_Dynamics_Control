% File LINZE.m
clear all
name = 'transp';
tfile = 'test';
%name = input('Enter Name of State Eqns. File : ‘,’s');
%tfile= input('Enter Name of Trim File : ‘,’s');
tmp= dlmread(tfile,',');
n=tmp(1); m=tmp(2); x=tmp(3:n+2);
u=tmp(n+3:m+n+2); tol=1e-6; time=0.;
mm = 4;
%mm= input('Number of control inputs to be used ? : ');
dx=0.1*x;
for i=1:n % Set Perturbations
if dx(i)==0.0
dx(i)=0.1;
end
end
last=zeros(n,1); a=zeros(n,n);
for j=1:n
xt=x;
for i=1:10
xt(j)=x(j)+dx(j);
xd1= feval (name,time,xt,u);
xt(j)=x(j)-dx(j);
xd2= feval (name,time,xt,u);
a(:,j)= (xd1-xd2)/(2*dx(j));
if max( abs(a(:,j)-last)./abs( a(:,j) + 1e-12 ) )<tol
break
end
dx(j)= 0.5*dx(j);
last = a(:,j);
end
%column=j
iteration=i;
disp(i)
if iteration==10
disp('not converged on A, column',num2str(j))
end
end
dlmwrite('A.dat',a);


