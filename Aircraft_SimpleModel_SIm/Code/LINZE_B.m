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

du=0.1*u;
for i=1:m % Set Perturbations
    if (i ~= m)
        if du(i)==0.0
        du(i)=0.1;
        end
    end
end
last=zeros(n,1); b=zeros(n,m);
for j=1:m
xt=x;
ut = u;
for i=1:10
xt(j)=x(j)+dx(j);
ut(j)=u(j)+du(j);
xd1= feval (name,time,xt,ut);
ut(j)=u(j)-du(j);
xd2= feval (name,time,xt,ut);


b(:,j)= (xd1-xd2)/(2*du(j));
if max( abs(b(:,j)-last)./abs( b(:,j) + 1e-12 ) )<tol
break
end
dx(j)= 0.5*dx(j);
du(j)= 0.5*du(j);
last = b(:,j);
end
%column=j
iteration=i;
if iteration==10
%disp('not converged on A, column',num2str(j))
end
end



