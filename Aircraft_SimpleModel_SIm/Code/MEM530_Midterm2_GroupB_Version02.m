%% Midterm 2 - Group B

%% Part a
clc; close all; clear all;
set(0,'defaultfigurecolor',[1 1 1])
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
title('Short-Period')
xlabel('Real'); ylabel('Imaginary')
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
title('Short-Period - No V')
xlabel('Real'); ylabel('Imaginary')
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
title('Phugoid Response')
xlabel('Real'); ylabel('Imaginary')
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
title('Phugoid Response - No V')
xlabel('Real'); ylabel('Imaginary')
grid on
hold on
plot([0 real(ph_noV(1))],[0 imag(ph_noV(1))])
hold on
plot([0 real(ph_noV(2))],[0 imag(ph_noV(2))])
hold on
plot([0 real(ph_noV(3))],[0 imag(ph_noV(3))])



%% Part c
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
initial(sys,x0,1000);
title('Velocity (ft/s) 10% Perturbed Phugoid')
subplot(2,1,2);
initial(sys,x0,20);
title('Velocity (ft/s) 10% Perturbed Short Period')
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
sys = canon(sys, 'modal')
figure(6)
subplot(2,1,1);
initial(sys,x0,1000)
title('Alpha (Radians) 10% Perturbed Phugoid')
subplot(2,1,2);

initial(sys,x0,20);
title('Alpha (Radians) 10% Perturbed Short Period')
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
figure(7)
subplot(2,1,1);
initial(sys,x0,1000)
title('Pitch (Radians) 10% Perturbed Phugoid')
subplot(2,1,2);
sys = canon(sys, 'modal')
initial(sys,x0,20);
title('Pitch (Radians) 10% Perturbed Short Period')
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
sys = canon(sys, 'modal')
figure(8)
subplot(2,1,1);
initial(sys,x0,1000)
title('Pitch Rate (Radians/Second) 10% Perturbed Phugoid')
subplot(2,1,2);
initial(sys,x0,20);
title('Pitch Rate (Radians/Second) 10% Perturbed Short Period')




%% Part d
%Plot time response using eigenvecotors
%Velocity
C = [1,0,0,0];
sys = ss(A,B,C,D);
figure(9)
subplot(2,2,1);
initial(sys,real(V(:,4)),1000)
title('Velocity Eigenvector Phugoid')
subplot(2,2,2);
initial(sys,real(V(:,2)),5)
title('Velocity  Eigenvector Short')
subplot(2,2,[3,4])
initial(sys,real(V(:,1))+real(V(:,3)),100)
title('Velocity Eigenvector Total')
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
figure(10)
subplot(2,2,1);
initial(sys,real(V(:,3)),1000)
title('Alpha Eigenvector Phugoid')
subplot(2,2,2);
initial(sys,real(V(:,1)),5)
title('Alpha Eigenvector Short')
subplot(2,2,[3,4])
initial(sys,real(V(:,4))+real(V(:,2)),100)
title('Alpha Eigenvector Total')
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
figure(11)
subplot(2,2,1);
initial(sys,real(V(:,4)),1000)
title('Pitch Eigenvector Phugoid')
subplot(2,2,2);
initial(sys,real(V(:,2)),5)
title('Pitch Eigenvector Short')
subplot(2,2,[3,4])
initial(sys,real(V(:,1))+real(V(:,3)),100)
title('Pitch Eigenvector Total')
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
figure(12)
subplot(2,2,1);
initial(sys,real(V(:,3)),1000)
title('Pitch Rate Eigenvector Phugoid')
subplot(2,2,2);
initial(sys,real(V(:,1)),5)
title('Pitch Rate Eigenvector Short')
subplot(2,2,[3,4])
initial(sys,real(V(:,2))+real(V(:,4)),100)
title('Pitch Rate Eigenvector Total')



%% Part e
%%
% Transfer function v/throttle
C = [1 0 0 0];
[NUM, DEN] = ss2tf(A,B,C,D,1);
sys_tf = tf(NUM, DEN)
figure(100)
subplot(2,1,1)
step(sys_tf,1000)
title('Velocity/Thrust Phugoid Step Response')
subplot(2,1,2)
step(sys_tf,5)
title('Velocity/Thrust Short Step Response')
figure(101)
rlocus(sys_tf)
title('Velocity/Thrust Root Locus')
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
%% 
%Transfer function a/elevator
C = [0 1 0 0];
[NUM, DEN] = ss2tf(A,B,C,D,2);
sys_tf = tf(NUM, DEN)
figure(102)
subplot(2,1,1)
step(sys_tf,1000)
title('Alpha/Elevator Phugoid Step Response')
subplot(2,1,2)
step(sys_tf,5)
title('Alpha/Elevator Short Step Response')
figure(103)
rlocus(sys_tf)
title('Alpha/Elevator Root Locus')
num = roots(NUM);
den = roots(DEN);
count = 0;
for i = 1:length(num)/2
    count2 = 0;
    for j = 1:length(den)/2
    if abs(abs(num(i+count)*num(i+count+1))-abs(den(j+count2)*den(j+count2+1))) < .1
        disp('There is a zero-pole cancelation!')
        Num = num(i+count)
        Den = den(j+count2)        
    end
    count2 = count2 + 1;
    end
    count = count + 1;
end

approx_sys_a = tf([0 0 -.04403], [1 2903/1250 22219369/5000000])
approx_sys_v = tf([0 8.2 -.18614], [1 49/5000 309313/50000000])
figure(104)
step(approx_sys_v,1000)
title('Velocity/Thrust Step Response Approx.')
figure(105)
step(approx_sys_a,10)
title('Alpha/Elevator Step Response Approx.')



%% Part f
% Generate the time-response plots of all states for (i) +0.1 units throttle with ‘zero’ elevator, and
% (ii) +2 deg. Elevator with ‘zero’ thrust. (These ‘inputs’ are perturbations from their trim
% values)
% Only two control surfaces for our model
trim_unp_input = [];
trim_unp_input = trim_u(:,1:2);
% Set time vectors for simulation
ph_t = 0:1:1000;
s_t = 0:1:1000;
% Add +0.1 Units to throttle, set elevator to zero
trim_p_throttle = [];
trim_p_throttle(1) = trim_unp_input(1) + 0.1;
trim_p_throttle(2) = 0.0;
% Add +2.2 Units to elevator, set throttle to zero
trim_p_elevator = [];
trim_p_elevator(1) = 0.0;
trim_p_elevator(2) = trim_unp_input(2) + 2.0;

%%
% Velocity Plots

% 'Velocity Response to Unpertubated Throttle and Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
C = [1 0 0 0];
sys = ss(A,B,C,D);
% 'Velocity Response to Pertubated Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_throttle(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);

figure(701)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Velocity Response to Pertubated Throttle and Zero Elevator','Throttle Value: ',num2str(trim_p_throttle(1))})

% 'Velocity Response to Unperturbed Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Velocity Response to Unperturbed Throttle and Zero Elevator','Throttle Value: ', num2str(trim_unp_input(1))})

% 'Velocity Response to Pertubated Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_p_elevator(2);

figure(702)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Velocity Response to Pertubated Elevator and Zero Throttle','Elevator Value: ',num2str(trim_p_elevator(2))})

% 'Velocity Response to Unperturbed Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Velocity Response to Unperturbed Elevator and Zero Throttle','Elevator Value: ', num2str(trim_unp_input(2))})

%%
% Alpha Plots

% 'Alpha Response to Unpertubated Throttle and Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
C = [0 1 0 0];
sys = ss(A,B,C,D);
% 'Alpha Response to Pertubated Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_throttle(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);

figure(704)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Alpha Response to Pertubated Throttle and Zero Elevator','Throttle Value: ',num2str(trim_p_throttle(1))})

% 'Alpha Response to Unperturbed Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Alpha Response to Unperturbed Throttle and Zero Elevator','Throttle Value: ', num2str(trim_unp_input(1))})

% 'Alpha Response to Pertubated Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_p_elevator(2);

figure(705)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Alpha Response to Pertubated Elevator and Zero Throttle','Elevator Value: ',num2str(trim_p_elevator(2))})

% 'Alpha Response to Unperturbed Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Alpha Response to Unperturbed Elevator and Zero Throttle','Elevator Value: ', num2str(trim_unp_input(2))})
%%
% Pitch Plots

% 'Pitch Response to Unpertubated Throttle and Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
C = [0 0 1 0];
sys = ss(A,B,C,D);
% 'Pitch Response to Pertubated Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_throttle(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);

figure(707)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Pitch Response to Pertubated Throttle and Zero Elevator','Throttle Value: ',num2str(trim_p_throttle(1))})

% 'Pitch Response to Unperturbed Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Pitch Response to Unperturbed Throttle and Zero Elevator','Throttle Value: ', num2str(trim_unp_input(1))})

% 'Pitch Response to Pertubated Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_p_elevator(2);

figure(708)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Pitch Response to Pertubated Elevator and Zero Throttle','Elevator Value: ',num2str(trim_p_elevator(2))})

% 'Pitch Response to Unperturbed Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Pitch Response to Unperturbed Elevator and Zero Throttle','Elevator Value: ', num2str(trim_unp_input(2))})

%%
% Pitch Rate Plots

% 'Pitch Rate Response to Unpertubated Throttle and Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
C = [0 1 0 0];
sys = ss(A,B,C,D);
% 'Pitch Rate Response to Pertubated Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_throttle(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);

figure(709)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Pitch Rate Response to Pertubated Throttle and Zero Elevator','Throttle Value: ',num2str(trim_p_throttle(1))})

% 'Pitch Rate Response to Unperturbed Throttle and Zero Elevator'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_unp_input(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Pitch Rate Response to Unperturbed Throttle and Zero Elevator','Throttle Value: ', num2str(trim_unp_input(1))})

% 'Pitch Rate Response to Pertubated Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_p_elevator(2);

figure(710)
subplot(2,2,[1,2])
lsim(sys,input_matrix,ph_t)
title({'Pitch Rate Response to Pertubated Elevator and Zero Throttle','Elevator Value: ',num2str(trim_p_elevator(2))})

% 'Pitch Rate Response to Unperturbed Elevator and Zero Throttle'
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_unp_input(2);
subplot(2,2,[3,4])
lsim(sys,input_matrix,ph_t)
title({'Pitch Rate Response to Unperturbed Elevator and Zero Throttle','Elevator Value: ', num2str(trim_unp_input(2))})

%% Part g
%% Part C, 10%
%Plot time response using eigenvecotors
%Velocity
x0 = [.1*x(1) ; .1*x(2) ; .1*x(3); .1*x(4)];
C = [1,0,0,0];
sys = ss(A,B,C,D);
[Vfc, tfv, xfv] = initial(sys,x0,1000);
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
[Afc, tfa, xfa] = initial(sys,x0,1000);
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
[Thfc, tfth, xfth] = initial(sys,x0,1000);
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
[Qfc, tfq, xfq] = initial(sys,x0,1000);

Vftc = x(1) + Vfc;
Aftc = x(2) + Afc;
Thftc = x(3) + Thfc;
gamc = zeros(1,length(Vfc)); Xc = zeros(1,length(Vfc)); hc = zeros(1,length(Vfc));
for i = 2:length(Vfc)
    gamc(i) = Thftc(i) - Aftc(i);
    Xc(i) = Xc(i-1)+(Vftc(i)*cos(gamc(i)))*(tfv(i)-tfv(i-1));
    hc(i) = hc(i-1)+(Vftc(i)*sin(gamc(i)))*(tfv(i)-tfv(i-1));
end
figure(1000)
plot(Xc,hc)
title('Flight Path of Vehicle (Inertial Frame) with 10% Perturbation In All States')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')


gamc = zeros(1,length(Vfc)); Xc = zeros(1,length(Vfc)); hc = zeros(1,length(Vfc));
for i = 2:length(Vfc)
    gamc(i) = Thftc(i) - Aftc(i);
    Xc(i) = Xc(i-1)+(Vftc(i)*cos(gamc(i))-x(1))*(tfv(i)-tfv(i-1));
    hc(i) = hc(i-1)+(Vftc(i)*sin(gamc(i)))*(tfv(i)-tfv(i-1));
end
figure(1001)
plot(Xc,hc)
title('Flight Path of Vehicle (Adjacent Vehicle Frame) with 10% Perturbation In All States')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')


%% Part D, Eigenvector
%Velocity
C = [1,0,0,0];
sys = ss(A,B,C,D);
[Vfd, tfv, xfv] = initial(sys,real(V(:,4)),1000);
%Alpha
C = [0 1 0 0];
sys = ss(A,B,C,D);
[Afd, tfa, xfa] = initial(sys,real(V(:,3)),1000);
%Pitch
C = [0 0 1 0];
sys = ss(A,B,C,D);
[Thfd, tfth, xfth] = initial(sys,real(V(:,4)),1000);
%Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
[Qfd, tfq, xfq] = initial(sys,real(V(:,3)),1000);

Vftd = x(1) + Vfd;
Aftd = x(2) + Afd;
Thftd = x(3) + Thfd;
gamd = zeros(1,length(Vfd)); Xd = zeros(1,length(Vfd)); hd = zeros(1,length(Vfd));
for i = 2:length(Vfd)
    gamd(i) = Thftd(i) - Aftd(i);
    Xd(i) = Xd(i-1)+(Vftd(i)*cos(gamd(i)))*(tfv(i)-tfv(i-1));
    hd(i) = hd(i-1)+(Vftd(i)*sin(gamd(i)))*(tfv(i)-tfv(i-1));
end
figure(1002)
plot(Xd,hd)
title('Flight Path of Vehicle (Inertial Frame) with Eigenvector Initial Conditions')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')


gamd = zeros(1,length(Vfd)); Xd = zeros(1,length(Vfd)); hd = zeros(1,length(Vfd));
for i = 2:length(Vfd)
    gamd(i) = Thftd(i) - Aftd(i);
    Xd(i) = Xd(i-1)+(Vftd(i)*cos(gamd(i))-x(1))*(tfv(i)-tfv(i-1));
    hd(i) = hd(i-1)+(Vftd(i)*sin(gamd(i)))*(tfv(i)-tfv(i-1));
end
figure(1003)
plot(Xd,hd)
title('Flight Path of Vehicle (Adjacent Vehicle Frame) with Eigenvector Initial Conditions')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')


%% f
%
% Velocity

% Input vector must be same length as the time vector
C = [1 0 0 0];
sys = ss(A,B,C,D);
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_throttle(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);
[Vff1, tff, xff] = lsim(sys,input_matrix,ph_t);

% Alpha
% Input vector must be same length as the time vector
C = [0 1 0 0];
sys = ss(A,B,C,D);
input_matrix = [];
input_matrix(1:length(s_t), 1) = trim_p_throttle(1);
input_matrix(1:length(s_t), 2) = trim_p_throttle(2);
[Aff1, tff, xff] = lsim(sys,input_matrix,s_t);

% Pitch
% Input vector must be same length as the time vector
C = [0 0 1 0];
sys = ss(A,B,C,D);
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_throttle(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);
[Thff1, tff, xff] = lsim(sys,input_matrix,ph_t);

% Pitch Rate
C = [0 0 0 1];
sys = ss(A,B,C,D);
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_throttle(1);
input_matrix(1:length(ph_t), 2) = trim_p_throttle(2);
[Qff1, tff, xff] = lsim(sys,input_matrix,ph_t);

Vf1t = x(1) + Vff1;
Af1t = x(2) + Aff1;
Thf1t = x(3) + Thff1;
gamf = zeros(1,length(Vff1)); Xf = zeros(1,length(Vff1)); hf = zeros(1,length(Vff1));
for i = 2:length(Vff1)
    gamf(i) = Thf1t(i) - Af1t(i);
    Xf(i) = Xf(i-1)+(Vf1t(i)*cos(gamf(i)))*(tff(i)-tff(i-1));
    hf(i) = hf(i-1)+(Vf1t(i)*sin(gamf(i)))*(tff(i)-tff(i-1));
end
figure(1004)
plot(Xf,hf)
title('Flight Path of Vehicle (Inertial Frame) with +0.1 Units Throttle Zero Elevator')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')


gamf = zeros(1,length(Vff1)); Xf = zeros(1,length(Vff1)); hf = zeros(1,length(Vff1));
for i = 2:length(Vff1)
    gamf(i) = Thf1t(i) - Af1t(i);
    Xf(i) = Xf(i-1)+(Vf1t(i)*cos(gamf(i))-x(1))*(tff(i)-tff(i-1));
    hf(i) = hf(i-1)+(Vf1t(i)*sin(gamf(i)))*(tff(i)-tff(i-1));
end
figure(1005)
plot(Xf,hf)
title('Flight Path of Vehicle (Adjacent Vehicle Frame) with +0.1 Units Throttle Zero Elevator')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')

% 'Velocity Response to Pertubated Elevator and Zero Throttle'
C = [1 0 0 0];
sys = ss(A,B,C,D);
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_p_elevator(2);
[Vff1, tff, xff] = lsim(sys,input_matrix,ph_t);

% 'Alpha Response to Pertubated Elevator and Zero Throttle'
C = [0 1 0 0];
sys = ss(A,B,C,D);
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_p_elevator(2);
[Aff1, tff, xff] = lsim(sys,input_matrix,ph_t);
% 'Pitch Response to Pertubated Elevator and Zero Throttle'
C = [0 0 1 0];
sys = ss(A,B,C,D);
input_matrix = [];
input_matrix(1:length(ph_t), 1) = trim_p_elevator(1);
input_matrix(1:length(ph_t), 2) = trim_p_elevator(2);
[Thff1, tff, xff] = lsim(sys,input_matrix,ph_t);

Vf1t = x(1) + Vff1;
Af1t = x(2) + Aff1;
Thf1t = x(3) + Thff1;
gamf = zeros(1,length(Vff1)); Xf = zeros(1,length(Vff1)); hf = zeros(1,length(Vff1));
for i = 2:length(Vff1)
    gamf(i) = Thf1t(i) - Af1t(i);
    Xf(i) = Xf(i-1)+(Vf1t(i)*cos(gamf(i)))*(tff(i)-tff(i-1));
    hf(i) = hf(i-1)+(Vf1t(i)*sin(gamf(i)))*(tff(i)-tff(i-1));
end
figure(1006)
plot(Xf,hf)
title('Flight Path of Vehicle (Inertial Frame) with Zero Throttle and +2deg Elevator')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')


gamf = zeros(1,length(Vff1)); Xf = zeros(1,length(Vff1)); hf = zeros(1,length(Vff1));
for i = 2:length(Vff1)
    gamf(i) = Thf1t(i) - Af1t(i);
    Xf(i) = Xf(i-1)+(Vf1t(i)*cos(gamf(i))-x(1))*(tff(i)-tff(i-1));
    hf(i) = hf(i-1)+(Vf1t(i)*sin(gamf(i)))*(tff(i)-tff(i-1));
end
figure(1007)
plot(Xf,hf)
title('Flight Path of Vehicle (Adjacent Vehicle Frame) with Zero Throttle and +2deg Elevator')
xlabel('Distance(ft)')
ylabel('Altitude (ft)')







