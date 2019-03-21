%% Part a
clc; close all; clear all;
[a, b] = linz_eqn('transp', 'test');
A = a(1:4, 1:4);
B = b(1:4, 1:2);
C = [1 0 0 0];
D = zeros(1,2);
sys = ss(A,B,C,D);



% Extract state vector and input vector from trim file
trim_read = dlmread('test', ',');
n = trim_read(1); % Number of states
m = trim_read(2); % Number of inputs
trim_x = trim_read(3:n+2); % State vector
trim_u=trim_read(n+3:m+n+2); % Input vector


% % Transfer function v/throttle
% [NUM, DEN] = ss2tf(A,B,C,D,1);
% sys_tf = tf(NUM, DEN)
% damp(sys_tf)
% 
% % Transfer function a/elevator
% [NUM, DEN] = ss2tf(A,B,C,D,2);
% C = [0 1 0 0];
% sys_tf = tf(NUM, DEN)
% damp(sys_tf)

[V, De] = eig(A);
% Short
eigenvalue = [De(1,1);De(2,2)];
freq = [(imag(De(1,1)));(imag(De(1,1)))];
damping_ratio = [-real(De(1,1))/(sqrt(real(De(1,1))^2+imag(De(1,1))^2)); -real(De(2,2))/(sqrt(real(De(2,2))^2+imag(De(2,2))^2))];
time_period = [2*pi/(abs(imag(De(1,1))));2*pi/(abs(imag(De(2,2))))];
Short = table(eigenvalue, freq, damping_ratio, time_period)
% Phugoid
eigenvalue = [De(3,3);De(4,4)];
freq = [(imag(De(3,3)));(imag(De(3,3)))];
damping_ratio = [-real(De(3,3))/(sqrt(real(De(3,3))^2+imag(De(3,3))^2));  -real(De(4,4))/(sqrt(real(De(4,4))^2+imag(De(4,4))^2))];
time_period = [2*pi/(abs(imag(De(3,3))));2*pi/(abs(imag(De(4,4))))];
Phugoid = table(eigenvalue, freq, damping_ratio, time_period)

%% Part b
% Short Period
sp = V(:,1);
labels = ["Vt" "a" "theta" "q"];
figure(1)
plot(real(sp),imag(sp),'*r') % marker plot
%text(real(sp),imag(sp),labels,'VerticalAlignment','top','HorizontalAlignment','left');
axis([-1 1 -1 1])
grid on
hold on
plot([0 real(sp(1))],[0 imag(sp(1))])
hold on
plot([0 real(sp(2))],[0 imag(sp(2))])
hold on
plot([0 real(sp(3))],[0 imag(sp(3))])
hold on
plot([0 real(sp(4))],[0 imag(sp(4))])

% Short Period w/o Velocity
sp_noV = [V(2,1) V(3,1) V(4,1)];
labels_noV = ["a" "theta" "q"];
figure(2)
plot(real(sp_noV),imag(sp_noV),'*r') % marker plot
%text(real(sp_noV),imag(sp_noV),labels_noV,'VerticalAlignment','top','HorizontalAlignment','left');
axis([-0.3 0.3 -0.3 0.3])
grid on
hold on
plot([0 real(sp_noV(1))],[0 imag(sp_noV(1))])
hold on
plot([0 real(sp_noV(2))],[0 imag(sp_noV(2))])
hold on
plot([0 real(sp_noV(3))],[0 imag(sp_noV(3))])

% Phugoid
ph = V(:,3);
figure(3)
plot(real(ph),imag(ph),'*r') % marker plot
%text(real(ph),imag(ph),labels,'VerticalAlignment','top','HorizontalAlignment','left');
axis([-1 1 -1 1])
grid on
hold on
plot([0 real(ph(1))],[0 imag(ph(1))])
hold on
plot([0 real(ph(2))],[0 imag(ph(2))])
hold on
plot([0 real(ph(3))],[0 imag(ph(3))])
hold on
plot([0 real(ph(4))],[0 imag(ph(4))])

% Phugoid w/o Velocity
ph_noV = [V(2,3) V(3,3) V(4,3)];
figure(4)
plot(real(ph_noV),imag(ph_noV),'*r') % marker plot
%text(real(ph_noV),imag(ph_noV),labels_noV,'VerticalAlignment','top','HorizontalAlignment','left');
axis([-0.003 0.003 -0.003 0.003])
grid on
hold on
plot([0 real(ph_noV(1))],[0 imag(ph_noV(1))])
hold on
plot([0 real(ph_noV(2))],[0 imag(ph_noV(2))])
hold on
plot([0 real(ph_noV(3))],[0 imag(ph_noV(3))])
%% Part c
%Plot the Response @ 10% i/c
trim_read = dlmread('test', ',');
n = trim_read(1); % Number of states
x = trim_read(3:n+2); % State vector
x0 = [.1*x(1) ; .1*x(2) ; .1*x(3); .1*x(4)];
%Velocity
C = [1 0 0 0];
sys = ss(A,B,C,D);
figure(5)
subplot(2,1,1);
initial(sys,x0);
title('V 10% Phugoid')
subplot(2,1,2);
% Converted to modal form for short period velocity
sys = canon(sys, 'modal');
initial(sys,x0,20);
title('V 10% S')
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
figure(6)
subplot(2,1,1);
initial(sys,x0)
title('A 10%')
subplot(2,1,2);
sys = canon(sys, 'modal');
initial(sys,x0,20);
title('A 10% S')
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
figure(7)
subplot(2,1,1);
initial(sys,x0)
title('Th 10%')
subplot(2,1,2);
sys = canon(sys, 'modal');
initial(sys,x0,20);
title('Th 10% S')
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
figure(8)
subplot(2,1,1);
initial(sys,x0)
title('Q 10%')
subplot(2,1,2);
sys = canon(sys, 'modal');
initial(sys,x0,20);
title('Q 10% S')
%% Part d
%Plot time response using eigenvecotors
%Velocity
C = [1,0,0,0];
sys = ss(A,B,C,D);
figure(9)
subplot(2,2,1);
initial(sys,real(V(:,4)))
title('V Eigenvectors')
subplot(2,2,2);
initial(sys,real(V(:,2)))
title('V Eigenvectors S')
subplot(2,2,[3,4])
initial(sys,real(V(:,1))+real(V(:,3)),100)
title('V Eigenvectors Total')
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
figure(10)
subplot(2,2,1);
initial(sys,real(V(:,3)))
title('A Eigenvectors')
subplot(2,2,2);
initial(sys,real(V(:,1)))
title('A Eigenvectors S')
subplot(2,2,[3,4])
initial(sys,real(V(:,4))+real(V(:,2)),100)
title('A Eigenvectors Total')
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
figure(11)
subplot(2,2,1);
initial(sys,real(V(:,4)))
title('Th Eigenvectors')
subplot(2,2,2);
initial(sys,real(V(:,2)))
title('Th Eigenvectors S')
subplot(2,2,[3,4])
initial(sys,real(V(:,1))+real(V(:,3)),100)
title('Th Eigenvectors Total')
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
figure(12)
subplot(2,2,1);
initial(sys,real(V(:,3)))
title('Q Eigenvectors')
subplot(2,2,2);
initial(sys,real(V(:,1)))
title('Q Eigenvectors S')
subplot(2,2,[3,4])
initial(sys,real(V(:,2))+real(V(:,4)),100)
title('Q Eigenvectors Total')
%% Part e
% Transfer function v/throttle
C = [1 0 0 0];
[NUM, DEN] = ss2tf(A,B,C,D,1);
sys_tf = tf(NUM, DEN)
figure(100)
subplot(2,1,1)
step(sys_tf)
title('V/Elevator')
subplot(2,1,2)
step(sys_tf,20)
title('V/Elevator S')
num = roots(NUM);
den = roots(DEN);
count = 0;
for i = 1:length(num)/2
    if abs(num(i+count)*num(i+count+1))-abs(den(i+count)*den(i+count+1)) < .1
        disp('There is a zero-pole cancelation!')
        Num = num(i+count)
        Den = den(i+count)
    end
    count = count + 1;
end
% Transfer function a/elevator
C = [0 1 0 0];
[NUM, DEN] = ss2tf(A,B,C,D,2);
sys_tf = tf(NUM, DEN)
figure(101)
subplot(2,1,1)
step(sys_tf)
title('A/Elevator')
subplot(2,1,2)
step(sys_tf,20)
title('A/Elevator S')
num = roots(NUM);
den = roots(DEN);
count = 0;
for i = 1:length(num)/2
    %FIX!!!
    if abs(abs(num(i+count)*num(i+count+1))-abs(den(i+count)*den(i+count+1))) < .1
        disp('There is a zero-pole cancelation!')
        Num = num(i+count)
        Den = den(i+count)
    end
    count = count + 1;
end

%% Part f
% Generate the time-response plots of all states for (i) +0.1 units throttle with ‘zero’ elevator, and
% (ii) +2 deg. Elevator with ‘zero’ thrust. (These ‘inputs’ are perturbations from their trim
% values)
% Only two control surfaces for our model
input_vector = trim_u(:,1:2);
% Add +0.1 Units to throttle, set elevator to zero
input_vector(1) = input_vector(1) + 0.1;
input_vector(2) = 0.0;
% Set time vector for simulation
t = 0:1:2500;
% Input vector must be same length as the time vector
input_matrix(1:length(t), 1) = input_vector(1);
input_matrix(1:length(t), 2) = input_vector(2);
% Get initial conditions
trim_x = trim_x(1:4);

%Velocity
C = [1 0 0 0];
sys = ss(A,B,C,D);
figure(200);
lsim(sys,input_matrix,t, trim_x')
title({'Velocity Response to Pertubated Throttle and Zero Elevator', num2str(input_vector(1))})

%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
figure(201);
lsim(sys,input_matrix,t, trim_x')
title({'Alpha Response to Pertubated Throttle and Zero Elevator', num2str(input_vector(1))})

%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
figure(202)
lsim(sys,input_matrix,t, trim_x')
title({'Pitch Response to Pertubated Throttle and Zero Elevator', num2str(input_vector(1))})

%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
figure(203);
lsim(sys,input_matrix,t, trim_x')
title({'Pitch Rate Response to Pertubated Throttle and Zero Elevator', num2str(input_vector(1))})
%%
% Only two control surfaces for our model
input_vector = trim_u(:,1:2);
% Add +2.2 Units to elevator, set throttle to zero
input_vector(1) = 0.0;
input_vector(2) = input_vector(2) + 2.0;
% Set time vector for simulation
t = 0:1:2500;
% Input vector must be same length as the time vector
input_matrix(1:length(t), 1) = input_vector(1);
input_matrix(1:length(t), 2) = input_vector(2);
% Get initial conditions
trim_x = trim_x(1:4);

%Velocity
C = [1 0 0 0];
sys = ss(A,B,C,D);
figure(204);
lsim(sys,input_matrix,t, trim_x')
title({'Velocity Response to Pertubated Elevator and Zero Throttle', num2str(input_vector(2))})

%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
figure(205);
lsim(sys,input_matrix,t, trim_x')
title({'Alpha Response to Pertubated Elevator and Zero Throttle', num2str(input_vector(2))})

%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
figure(206)
lsim(sys,input_matrix,t, trim_x')
title({'Pitch Response to Pertubated Elevator and Zero Throttle', num2str(input_vector(2))})

%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
figure(207);
lsim(sys,input_matrix,t, trim_x')
title({'Pitch Rate Response to Pertubated Elevator and Zero Throttle', num2str(input_vector(2))})

%% Part g
%% i
% Part C, 10%
%Plot time response using eigenvecotors
%Velocity
x0 = [.1*x(1) ; .1*x(2) ; .1*x(3); .1*x(4)];
C = [1,0,0,0];
sys = ss(A,B,C,D);
[Vf, tfv, xfv] = initial(sys,x0,3000);
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
[Af, tfa, xfa] = initial(sys,x0,3000);
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
[Thf, tfth, xfth] = initial(sys,x0,3000);
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
[Qf, tfq, xfq] = initial(sys,x0,3000);

Vftc = x(1) + Vf;
Aftc = x(2) + Af;
Thftc = x(3) + Thf;
gam = zeros(1,length(Vf)); X = zeros(1,length(Vf)); h = zeros(1,length(Vf));
for i = 2:length(Vf)
    gam(i) = Thftc(i) - Aftc(i);
    X(i) = X(i-1)+(Vftc(i)*cos(gam(i)))*(tfv(i)-tfv(i-1));
    h(i) = h(i-1)+(Vftc(i)*sin(gam(i)))*(tfv(i)-tfv(i-1));
end
figure(1000)
plot(X,h)

gam = zeros(1,length(Vf)); X = zeros(1,length(Vf)); h = zeros(1,length(Vf));
for i = 2:length(Vf)
    gam(i) = Thftc(i) - Aftc(i);
    X(i) = X(i-1)+(Vftc(i)*cos(gam(i))-x(1))*(tfv(i)-tfv(i-1));
    h(i) = h(i-1)+(Vftc(i)*sin(gam(i)))*(tfv(i)-tfv(i-1));
end
figure(1001)
plot(X,h)

% Part D, Eigenvector
%Velocity
C = [1,0,0,0];
sys = ss(A,B,C,D);
[Vf, tfv, xfv] = initial(sys,real(V(:,4)),3000);
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
[Af, tfa, xfa] = initial(sys,real(V(:,3)),3000);
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
[Thf, tfth, xfth] = initial(sys,real(V(:,4)),3000);
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
[Qf, tfq, xfq] = initial(sys,real(V(:,3)),3000);

Vft = x(1) + Vf;
Aft = x(2) + Af;
Thft = x(3) + Thf;
gam = zeros(1,length(Vf)); X = zeros(1,length(Vf)); h = zeros(1,length(Vf));
for i = 2:length(Vf)
    gam(i) = Thft(i) - Aft(i);
    X(i) = X(i-1)+(Vft(i)*cos(gam(i)))*(tfv(i)-tfv(i-1));
    h(i) = h(i-1)+(Vft(i)*sin(gam(i)))*(tfv(i)-tfv(i-1));
end
figure(1002)
plot(X,h)

gam = zeros(1,length(Vf)); X = zeros(1,length(Vf)); h = zeros(1,length(Vf));
for i = 2:length(Vf)
    gam(i) = Thft(i) - Aft(i);
    X(i) = X(i-1)+(Vft(i)*cos(gam(i))-x(1))*(tfv(i)-tfv(i-1));
    h(i) = h(i-1)+(Vft(i)*sin(gam(i)))*(tfv(i)-tfv(i-1));
end
figure(1003)
plot(X,h)
