%% Problem 1 Part A: Determine Poles and Zeros of Transfer Function
% Generate transfer function in matlab
num = [-0.01785 -1.38873 -0.0088536 -0.0079611]
denom = [1 0.81066 1.32005596 0.01038106 0.0069483]
tfsys = tf(num, denom)
%% 
% Compute the poles of the transfer function model
poles = pole(tfsys)
%%
% Compute the zeros of the transfer function model
zeros = zero(tfsys)
%% Problem 1 Part B: Is the system stable?
% Looking at the poles and zeros we can generate
% pole-zero plot. For BIBO stability the poles must be in the left hand
% plane
pzmap(tfsys)
grid on
%%
% We can see here that our poles are in the left hand plane, so this
% system is BIBO stable. BIBO stability also states that the roots of our
% transfer function must have negative real parts. We can observe:
real(poles)
%%
% That our values are indeed negative real, as well as repeated roots.
% Finally, we can plot an unit impulse response of our system to see if it
% goes to zero as time approaches infinity
impulse(tfsys)
%%
% Our system reaches zero as time moves to infinity, so its BIBO stable in
% all cases
%%
% If we looked into internal stability, we see that with our repeated roots
% of our transfer function our system is not marginally stable. Since the real roots have a negative
% part however, we are asymptotically stable. 
%%
% What are our natural frequencies and damping ratio?
[W, zeta] = damp(tfsys)
%% Problem 1 Part D: Convert to State Space Form
[A, B, C, D] = tf2ss(num, denom)
sys = ss(A, B, C, D)
%% Problem 1 Part E: Plot step response
step(sys)
%% Problem 2 Part A: Determine Poles and Zeros of Transfer Function
% Generate transfer function in matlab
num2= [-0.01782 -1.386396]
denom2 = [1 0.805 1.325]
tfsys2 = tf(num2, denom2)
%% 
% Compute the poles of the transfer function model
poles2 = pole(tfsys2)
%%
% Compute the zeros of the transfer function model
zeros2 = zero(tfsys2)
%% Problem 2 Part B: Is the system stable?
% Looking at the roots of the denominator (the poles) we can generate
% pole-zero plot. For BIBO stability the poles must be in the left hand
% plane
pzmap(tfsys2)
grid on
%%
% We can see here that our poles are in the left hand plane, so this
% system is BIBO stable. BIBO stability also states that the roots of our
% transfer function must have negative real parts. We can observe:
real(poles2)
%%
% That our values are indeed negative real, as well as repeated roots.
% Finally, we can plot an unit impulse response of our system to see if it
% goes to zero as time approaches infinity
impulse(tfsys2)
%%
% Our system reaches zero as time moves to infinity, so its BIBO stable in
% all cases
%%
% If we looked into internal stability, we see that with our repeated roots
% of our transfer function our system is not marginally stable. Since the real roots have a negative
% part however, we are asymptotically stable. 
%%
% What are our natural frequencies and damping ratio?
[W, zeta] = damp(tfsys)
%% Problem 2 Part D: Convert to State Space Form
[A2, B2, C2, D2] = tf2ss(num2, denom2)
sys2 = ss(A2, B2, C2, D2)
%% Problem 2 Part E: Plot step response
step(sys2)
%% Problem 3 Compare both systems
% We can compare the responses between both systems
step(sys)
hold on;
step(sys2)
legend('System 1', 'System 2')
%% Problem 4 Part A: Moment of Inertia of Body without Tip Masses
% J = (1/12)*M*(h^2 + w^2) for xaxis
%% Setup Variables
m_body = 10 %kg
m_tips = 2 %kg
d_body = {10,10,30} %cm
d_tips = {5, 5, 10} %cm
y = 25 %cm
J_body = (1/12)*m_body*( (d_body{1}^2) + (d_body{2}^2) ) %In kg*cm^2
%% Problem 4 Part B: Moment of Inertia for Tips
J_tips = (1/12)*m_tips*( (d_tips{1}^2) + (d_tips{2}^2) ) %In kg*cm^2
%% Problem 4 Part C: Moment of Inertia for Entire System
J_total = (J_body) + ((m_tips*(y^2))*2)