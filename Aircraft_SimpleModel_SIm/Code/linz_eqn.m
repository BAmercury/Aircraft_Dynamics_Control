function [a,b] = linz_eqn(model, trim_values)
%   model: string file name of the model (example: 'transp')
%   trim_values: string file name of the trim value output file (example:
%   'test')
%   Typical application: 
%   [a,b] = linz_eqn('transp', 'test')

% Extract state vector and input vector from trim file
trim_read = dlmread(trim_values, ',');
n = trim_read(1); % Number of states
m = trim_read(2); % Number of inputs
x = trim_read(3:n+2); % State vector
u=trim_read(n+3:m+n+2); % Input vector
tol=1e-6; time=0.; % Set tolerance and time

%%
% Set Pertubations for x (state vector), default is 10%
level =0.1;
dx = level*x;
for i=1:n
    if dx(i) == 0.0
        dx(i) = level;
    end
end

% Set Pertubations for u (input vector)
du=level*u;
for i=1:m % Set Perturbations
    if (i ~= m) % Make sure it skips over landing gear input value
        if du(i)==0.0
        du(i)=level;
        end
    end
end
%%

% Create empty a and b matrices
lasta = zeros(n,1);
a = zeros(n,n);

lastb = zeros(n,1);
b = zeros(n,m);

%%
% Iterate and generate columns of a
for j=1:n
    xt=x;
    for i=1:10 % Step size 10
        xt(j)=x(j)+dx(j);
        xd1= feval (model,time,xt,u);
        xt(j)=x(j)-dx(j);
        xd2= feval (model,time,xt,u);
        a(:,j)= (xd1-xd2)/(2*dx(j));
        if max( abs(a(:,j)-lasta)./abs( a(:,j) + 1e-12 ) )<tol
            break
        end
        dx(j)= 0.5*dx(j);
        lasta = a(:,j);
    end
end
    
% now iterate and generate the columns of b
for j=1:m
    xt=x;
    ut = u;
    for i=1:10
        xt(j)=x(j)+dx(j);
        ut(j)=u(j)+du(j);
        xd1= feval (model,time,xt,ut);
        ut(j)=u(j)-du(j);
        xd2= feval (model,time,xt,ut);


        b(:,j)= (xd1-xd2)/(2*du(j));
        if max( abs(b(:,j)-lastb)./abs( b(:,j) + 1e-12 ) )<tol
            break
        end
        dx(j)= 0.5*dx(j);
        du(j)= 0.5*du(j);
        lastb = b(:,j);
    end

end



end