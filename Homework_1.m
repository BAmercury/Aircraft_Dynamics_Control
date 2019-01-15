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
% Looking at the roots of the denominator (the poles) we can generate
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
% part however, we are asympototically stable. 

%%
% What are our natural frequencies and damping ratio?
[W, zeta] = damp(tfsys)

%% Problem 1 Part C: Are there any pole-zero cancellations?
% If we take our transfer function back to state space, and observe
% eigenvalues of our A matrix. We can see if any of the poles or zeros were
% cancelled
[A, B, C, D] = tf2ss(num, denom)
%%
% Take our realization and put it into state space form
sys = ss(A, B, C, D)
%%
% Find poles and zeros
ss_poles = pole(sys)
ss_zeros = zero(sys)
%%
% Check eigenvalues of A as well to see if it matches with poles of our
% transfer function
eig(A)
%%
% As we can see, there is no pole-zero cancellation

%% Problem 1 Part D: Convert to State Space Form
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

% We can compare this with our original system
hold on;
impulse(tfsys)
%%
% Our system reaches zero as time moves to infinity, so its BIBO stable in
% all cases
%%
% If we looked into internal stability, we see that with our repeated roots
% of our transfer function our system is not marginally stable. Since the real roots have a negative
% part however, we are asympototically stable. 

%%
% What are our natural frequencies and damping ratio?
[W, zeta] = damp(tfsys)

%% Problem 2 Part C: Are there any pole-zero cancellations?
% If we take our transfer function back to state space, and observe
% eigenvalues of our A matrix. We can see if any of the poles or zeros were
% cancelled
[A2, B2, C2, D2] = tf2ss(num2, denom2)
%%
% Take our realization and put it into state space form
sys2 = ss(A2, B2, C2, D2)
%%
% Find poles and zeros
ss_poles = pole(sys2)
ss_zeros = zero(sys2)
%%
% Check eigenvalues of A as well to see if it matches with poles of our
% transfer function
eig(A2)
%%
% As we can see, there is no pole-zero cancellation
%% Problem 2 Part D: Convert to State Space Form
sys2 = ss(A2, B2, C2, D2)

%% Problem 2 Part E: Plot step response
step(sys2)
hold on;
step(sys)


