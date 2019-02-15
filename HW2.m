%% MEM 530 Homework 2
% Bhautik (Brian) Amin


%% Problem 1
clear
clc
% Brick properties
prop = [8;5;2]; % Length, Width, Height (X,Y,Z)
mass = 12; 
% Calcuate moment of inertia
J = zeros(3,1);
J(1) = (mass/12) * ((prop(3)^2) + (prop(2)^2)); % X axis
J(2) = (mass/12) * ((prop(3)^2) + (prop(1)^2)); % Y axis
J(3) = (mass/12) * ((prop(2)^2) + (prop(1)^2)); % Z axis


tspam = [0:0.001:100]; % Simulate for 100 seconds

% Set up initial conditions
parta_ic = 180/pi * [0.1; 0; 0.001]; % deg/s
partb_ic = 180/pi * [0.001; 0; 0.1];
partc_ic = 180/pi * [0; 0.1; 0.001];
%partc_ic = 180/pi * [0; 0; 0];

% ODE Solver Part A
[t,y] = ode45(@(t,y) euler_motion_ode(t,y,J), tspam, parta_ic);
figure(1);
plot(t,y(:,1),'b',t,y(:,2),'g',t,y(:,3),'r')
title({'Part A Initial Conditions','P= 5.7296,Q=0,R=0.0573 deg/s'})
xlabel('Time t')
ylabel('Angular Velocity (Deg/s)')
legend('P','Q','R')

% ODE Solver Part B
[t,y] = ode45(@(t,y) euler_motion_ode(t,y,J), tspam, partb_ic);
figure(2);
plot(t,y(:,1),'b',t,y(:,2),'g',t,y(:,3),'r')
title({'Part B Initial Conditions','P= 0.0573,Q=0,R=5.7296 deg/s'})
xlabel('Time t')
ylabel('Angular Velocity')
legend('P','Q','R')

% ODE Solver Part C
[t,y] = ode45(@(t,y) euler_motion_ode(t,y,J), tspam, partc_ic);
figure(3);
plot(t,y(:,1),'b',t,y(:,2),'g',t,y(:,3),'r')
title({'Part C Initial Conditions','P=0,Q=5.7296,R=0.0573 deg/s'})
xlabel('Time t')
ylabel('Angular Velocity')
legend('P','Q','R')
