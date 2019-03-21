%% MEM 530 Final - Team B

%% Step 1

%% Given Airplane Values

clear; clc; close all;
S = 5170; m = 5000; Jy = 4.1*10^6;
cbar = 17.5; VTe = 500; qbar = 297.125; 
Cma_dot = -6; Cma = -0.022; Cmq = -16; Cmde = -0.016;
alpha = 0.5798; Cd = 0.016 + .042*(.2 + .085*alpha)^2;
Cla = 0.085; Clde = 0;
Za = -qbar*S/m*(Cd + Cla); Za_dot = 0; Zq = 0; Zde = -qbar*S/m*Clde;
Ma_dot = (qbar*S*cbar/Jy)*(cbar/2/VTe)*Cma_dot;
Ma = (qbar*S*cbar/Jy)* Cma; Mq = (qbar*S*cbar/Jy)*(cbar/2/VTe)*Cmq;
Mde = (qbar*S*cbar/Jy)* Cmde;
% E*xdot = A*x + B*u , E = [VTe - Za_dot, 0; -Ma_dot, 1]

%% Formation of the Short Period Model

%A = [-0.2082 -128.6490  -32.1700 0; -0.0017   -1.2600 \0  1.0000; 0   0   0  1.0000; -0.0070   -3.1048  0 -1.0595];  %Ashort = [-1.2608, 1.00; -3.1046, -1.0595]; Ashort(:,3) = 0;Ashort(1,:) = 0;Ashort(3,:) = 0;
A = [-0.0118, 18.9953 , -32.1700,  0; -0.0003, -1.2608,  0, 1.0000;  0   0  0  1.0000;    0.0000   -3.1046         0   -1.0595];
Ashort = [A(2,2) A(2,4); A(4,2) A(4,4)];
%Bshort = [0.02; -0.0002];
B =[8.1996 0; -0.0002 0; 0 0; 0.0200   -0.0440];
Bshort = B(2:2:4,2);
C = [180/pi 0; 0 180/pi];
D = 0;
% Test
% Ashort = [-.019311 8.8157 -32.17 -.574999; -2.5389*10^-4 -1.0189 0 .90506; 0 0 0 1; 2.9465*10^-12 .82225 0 -1.0774];
% Bshort= [.1737; -2.1499*10^-3; 0; -.17555];
% C = [0 57.29578 0 0; 0 0 0 57.29578];
for j = 1:2
    Ga(j,:) = [Ashort(j,:) -Bshort(j,1) 0];
end
Ga(3,:) = [0 0 -20.2 0];
Ga(4,:) = [10 0 0 -10];
Gb = [0;0;20.2;0];
Gc = [C(1,:) 0 0; C(2,:) 0 0; 0 0 0 180/pi];
Gd = 0;



%% Step 2

%% Determine the AOA Gain

sys = ss(Ga,Gb,Gc(3,:),Gd);
figure(1)
rlocus(sys)
disp('========== SHORT-PERIOD MODEL FLIGHT WITH NO SAS ==========') 
[NUM,DEN] = ss2tf(Ga,Gb,Gc(3,:),Gd);
Gshort_total = tf(NUM,DEN)
root = (roots(DEN));
real = real(root(3));
imag = imag(root(3));
Freq1 = sqrt(real^2+imag^2)
Damp1 =  -real/(sqrt(real^2+imag^2))
zeta = 0.95; wn = 4.8; sgrid(zeta,wn) % [k,poles] = rlocfind(sys)
ka = 2.938;

%% Determine the Pitch Gain

aca = Ga-Gb*ka*Gc(3,:);
sysa = ss(aca,Gb,Gc(2,:),Gd);
figure(2)
rlocus(sysa);
disp('========== SHORT-PERIOD MODEL FLIGHT WITH Ka IMPLEMENTED ==========') 
[NUMa,DENa] = ss2tf(aca,Gb,Gc(2,:),Gd);
Gshort_a_total = tf(NUMa,DENa)
% Find freuency and damping ratio
root = (roots(DENa));
Freq2 = sqrt(.6013^2+3.1808^2)
Damp2 =  -real/(sqrt(.6013^2+3.1808^2))
zeta = 0.95;  wn = 4.8;  sgrid(zeta,wn)
kq = 2.7676;

%% Final Transfer Function

acq = aca-Gb*kq*Gc(2,:);
sysq = ss(acq,Gb,Gc(2,:),Gd);
figure(3)
rlocus(sysq);
disp('========== SHORT-PERIOD MODEL FLIGHT WITH Ka AND Kq IMPLEMENTED ==========') 
[NUMq,DENq] = ss2tf(acq,Gb,Gc(2,:),0);
Gshort_q_total = tf(NUMq,DENq)
root = (roots(DENq));
Freq3 = sqrt(4.5602^2+1.498^2)
Damp3 =  4.5602/(sqrt(4.5602^2+1.498^2))
zeta = 0.95;  wn = 4.8;  sgrid(zeta,wn)


%% Step 3

%% Ploting with Initial Response

trim_read = dlmread('test', ',');
n = trim_read(1); % Number of states
x = trim_read(3:n+2); % State vector
%a
sys = ss(acq,Gb,Gc(1,:),D);
x0 = [.1*x(2) ; .1*x(4); 0; 0];
figure(4)
initial(sys,x0);
ylabel('Amplitude (deg)')
title('Short Period SAS - Alpha 10% Perturbed')
%q
sys = ss(acq,Gb,Gc(2,:),D);
x0 = [.1*x(2) ; .1*x(4); 0; 0];
figure(5)
initial(sys,x0);
ylabel('Amplitude (deg/sec)')
title('Short Period SAS - Pitch Rate 10% Perturbed')
%  a elevator step response simulation
sys = ss(acq,Gb,Gc(1,:),D);
figure(100)
step(sys);
ylabel('Amplitude (deg)')
xlabel('Time (seconds')
title('Elevator Step Response - Alpha')
% q elevator step response simulation
sys = ss(acq,Gb,Gc(2,:),D);
figure(101)
step(sys)
ylabel('Amplitude (deg/sec)')
xlabel('Time (seconds)')
title('Elevator Step Response - Pitch Rate')



%% Step 4

%% Full Longitudinal State Space

A = [-0.0118, 18.9953 , -32.1700,  0; -0.0003, -1.2608,  0, 1.0000;  0   0  0  1.0000;    0.0000   -3.1046         0   -1.0595];
B =[8.1996 0; -0.0002 0; 0 0; 0.0200   -0.0440];
C = [0 180/pi 0 0; 0 0 0 180/pi];
D = 0; 
Ga = zeros(6,6);
for j = 1:4
    Ga(j,:) = [A(j,:) -B(j,2) 0];
end
Ga(5,:) = [0 0 0 0 -20.2 0];
Ga(6,:) = [0 10 0 0 0 -10];
Gb = [0; 0; 0; 0; 20.2; 0];
Gc = [C(1,:) 0 0; C(2,:) 0 0; 0 0 0 0 0 180/pi];
Gd = 0;

%% Determine the AOA Gain and Pitch Gain

sys = ss(Ga,Gb,Gc(3,:),Gd);
figure(6)
rlocus(sys)
[NUM,DEN] = ss2tf(Ga,Gb,Gc(3,:),Gd);
Glong_total = tf(NUM,DEN);
ka = 2.938;

aca = Ga-Gb*ka*Gc(3,:);
sysa = ss(aca,Gb,Gc(2,:),Gd);
figure(7)
rlocus(sysa);
[NUMa,DENa] = ss2tf(aca,Gb,Gc(2,:),Gd);
Gshort_a_total = tf(NUMa,DENa);
kq = 2.7676;

%% Final Transfer Function

acq = aca-Gb*kq*Gc(2,:);
sysq = ss(acq,Gb,Gc(2,:),Gd);
figure(8)
rlocus(sysq);
disp('========== FULL LONGITUDINAL FLIGHT WITH SAS IMPLEMENTED ==========') 
[NUMq,DENq] = ss2tf(acq,Gb,Gc(2,:),0);
Gshort_q_total = tf(NUMq,DENq)
root = (roots(DENq));
Freq4 = sqrt(4.5602^2+1.498^2)
Damp4 =  4.5602/(sqrt(4.5602^2+1.498^2))

%% Plot Initial Response

trim_read = dlmread('test', ',');
n = trim_read(1); % Number of states
x = trim_read(3:n+2); % State vector
x0 = [.1*x(1) ; .1*x(2) ; .1*x(3); .1*x(4); 0; 0];
%Velocity
C = [1 0 0 0 0 0];
sys = ss(acq,Gb,C,D);
figure(9)
subplot(2,1,1);
initial(sys,x0);
ylabel('Amplitude (ft/sec)')
title('Full Longitudinal SAS - Velocity 10% Perturbed')
subplot(2,1,2);
initial(sys,x0,.2);
ylabel('Amplitude (ft/sec)')
title('Zoomed In To Show Short Period')
%Alpha
C = [0 180/pi 0 0 0 0];
sys = ss(acq,Gb,C,D);
figure(10)
subplot(2,1,1);
initial(sys,x0);
ylabel('Amplitude (deg)')
title('Full Longitudinal SAS - Alpha 10% Perturbed')
subplot(2,1,2);
initial(sys,x0,.2);
ylabel('Amplitude (deg)')
title('Zoomed In To Show Short Period')
%Pitch
C = [ 0 0 180/pi 0 0 0];
sys = ss(acq,Gb,C,D);
figure(11)
subplot(2,1,1);
initial(sys,x0);
ylabel('Amplitude (deg)')
title('Full Longitudinal SAS - Pitch 10% Perturbed')
subplot(2,1,2);
initial(sys,x0,.2);
ylabel('Amplitude (deg)')
title('Zoomed In To Show Short Period')
%Pitch Rate
C = [0 0 0 180/pi 0 0];
sys = ss(acq,Gb,C,D);
figure(12)
subplot(2,1,1);
initial(sys,x0);
ylabel('Amplitude (deg/sec)')
title('Full Longitudinal SAS - Pitch Rate 10% Perturbed')
subplot(2,1,2);
initial(sys,x0,.2);
ylabel('Amplitude (deg/sec)')
title('Zoomed In To Show Short Period')
%%
%Velocity Elevator Step Response
C = [1 0 0 0 0 0];
sys = ss(acq,Gb,C,D);
figure(200)
step(sys)
ylabel('Amplitude (ft/sec)')
xlabel('Time (seconds)')
title({'Full Longitudinal SAS Elevator Step Response', 'Velocity'})

% Alpha Elevator Step Response
C = [0 180/pi 0 0 0 0];
sys = ss(acq,Gb,C,D);
figure(201)
step(sys)
ylabel('Amplitude (deg)')
xlabel('Time (seconds)')
title({'Full Longitudinal SAS Elevator Step Response', 'Alpha'})

% Pitch Elevator Step Response
C = [ 0 0 180/pi 0 0 0];
sys = ss(acq,Gb,C,D);
figure(202)
step(sys)
ylabel('Amplitude (deg)')
xlabel('Time (seconds)')
title({'Full Longitudinal SAS Elevator Step Response', 'Pitch'})


% Pitch Rate Elevator Step Response
C = [0 0 0 180/pi 0 0];
sys = ss(acq,Gb,C,D);
figure(203)
step(sys)
ylabel('Amplitude (deg/s)')
xlabel('Time (seconds)')
title({'Full Longitudinal SAS Elevator Step Response', 'Pitch Rate'})



%% Part 5

%% Formation of the Short Period Model

%% Short Period Model without SAS
A = [-0.0118, 18.9953 , -32.1700,  0; -0.0003, -1.2608,  0, 1.0000;  0   0  0  1.0000;    0.0000   -3.1046         0   -1.0595];
Ashort = [A(2,2) A(2,4); A(4,2) A(4,4)];
B =[8.1996 0; -0.0002 0; 0 0; 0.0200   -0.0440];
Bshort = B(2:2:4,2);
C = [180/pi 0; 0 180/pi];
D = 0;
for j = 1:2
    Ga(j,:) = [Ashort(j,:) -Bshort(j,1) 0];
end
Ga(3,:) = [0 0 -20.2 0];
Ga(4,:) = [10 0 0 -10];
Gb = [0;0;20.2;0];
Gc = [C(1,:) 0 0; C(2,:) 0 0; 0 0 0 180/pi];
Gd = 0;
%%
sys = ss(Ga,Gb,Gc(3,:),Gd);
disp('========== SHORT-PERIOD MODEL FLIGHT WITH NO SAS ==========') 
[NUM,DEN] = ss2tf(Ga,Gb,Gc(3,:),Gd);
Gshort_no_sas = tf(NUM,DEN)
%%
%% Ploting with Initial Response
trim_read = dlmread('test', ',');
n = trim_read(1); % Number of states
x = trim_read(3:n+2); % State vector
%a
sys = ss(Ga,Gb,Gc(1,:),Gd);
x0 = [.1*x(2) ; .1*x(4); 0; 0];
figure(400)
initial(sys,x0);
ylabel('Amplitude (deg)')
title({'Short Period SAS and No SAS - Alpha 10% Perturbed'})
hold on;
sys = ss(acq,Gb,Gc(1,:),D);
initial(sys,x0);
legend('No SAS', 'SAS')
%q
sys = ss(Ga,Gb,Gc(2,:),Gd);
x0 = [.1*x(2) ; .1*x(4); 0; 0];
figure(500)
initial(sys,x0);
ylabel('Amplitude (deg/sec)')
title({'Short Period SAS and No SAS - Pitch Rate 10% Perturbed'})
hold on;
sys = ss(acq,Gb,Gc(2,:),D);
initial(sys,x0);
legend('No SAS', 'SAS')



%%
% Full longitudal no sas vs sas
%% Full Longitudinal State Space
A = [-0.0118, 18.9953 , -32.1700,  0; -0.0003, -1.2608,  0, 1.0000;  0   0  0  1.0000;    0.0000   -3.1046         0   -1.0595];
B =[8.1996 0; -0.0002 0; 0 0; 0.0200   -0.0440];
C = [0 180/pi 0 0; 0 0 0 180/pi];
D = 0; 
Ga = zeros(6,6);
for j = 1:4
    Ga(j,:) = [A(j,:) -B(j,2) 0];
end
Ga(5,:) = [0 0 0 0 -20.2 0];
Ga(6,:) = [0 10 0 0 0 -10];
Gb = [0; 0; 0; 0; 20.2; 0];
Gc = [C(1,:) 0 0; C(2,:) 0 0; 0 0 0 0 0 180/pi];
Gd = 0;
%%
sys = ss(Ga,Gb,Gc(3,:),Gd);
ka = 2.938;
aca = Ga-Gb*ka*Gc(3,:);
kq = 2.7676;
acq = aca-Gb*kq*Gc(2,:);

%% Ploting with Initial Response
trim_read = dlmread('test', ',');
n = trim_read(1); % Number of states
x = trim_read(3:n+2); % State vector
x0 = [.1*x(1) ; .1*x(2) ; .1*x(3); .1*x(4); 0; 0];
%Velocity
figure(600);
C = [1 0 0 0 0 0];
sys = ss(Ga,Gb,C,D);
initial(sys,x0);
ylabel('Amplitude (ft/sec)')
title('Full Longitudinal No SAS vs SAS - Velocity 10% Perturbed')
hold on;
sys = ss(acq,Gb,C,D);
initial(sys,x0);
legend('No SAS', 'SAS')

%Alpha
C = [0 180/pi 0 0 0 0];
sys = ss(Ga,Gb,C,D);
figure(601)
initial(sys,x0);
ylabel('Amplitude (deg)')
title('Full Longitudinal No SAS vs SAS - Alpha 10% Perturbed')
hold on;
sys = ss(acq,Gb,C,D);
initial(sys,x0);
legend('No SAS', 'SAS')

%Pitch
C = [ 0 0 180/pi 0 0 0];
sys = ss(Ga,Gb,C,D);
figure(602)
initial(sys,x0);
ylabel('Amplitude (deg)')
title('Full Longitudinal No SAS vs SAS - Pitch 10% Perturbed')
hold on;
sys = ss(acq,Gb,C,D);
initial(sys,x0);
legend('No SAS', 'SAS')



%Pitch Rate
C = [0 0 0 180/pi 0 0];
sys = ss(Ga,Gb,C,D);
figure(603)
initial(sys,x0);
ylabel('Amplitude (deg/sec)')
title('Full Longitudinal No SAS vs SAS - Pitch Rate 10% Perturbed')
hold on;
sys = ss(acq,Gb,C,D);
initial(sys,x0);
legend('No SAS', 'SAS')