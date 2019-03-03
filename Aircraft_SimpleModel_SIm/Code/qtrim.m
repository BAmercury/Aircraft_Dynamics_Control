%% quick trim
global rho
global num_1 num_2 den_1 den_2 den_3
global cdcls cl0 cla
s = 2170; cbar = 17.5; mass = 5.0e+3; Iyy = 4.1E06; 
tstat = 6.0e4; dtdv = -38.0; ze = 2.0; cdcls = .042;
cla = .085; cma = -.022; cmde = -.016; % per degree
cmq = -16.0; cmadot = -6.0; cladot = 0.0; % per radian
cl0 = 0.2; cd0 = 0.016;
cm0 = 0.05; dcdg = 0.0; dcmg = 0.0;
rtod = 57.29578; gd = 32.17;

%        ze = 0.0;   %% assumption

%% define eqlbm flt
vt = input('velocity (ft/sec.) ? : ');
h = input('altitude (ft.) ? : ');
gam = input('climb angle (deg.) ? : ')/rtod;        %% gam in radians

%% arbitrarily set alpha to zero (level fuselage)
alpha = 0;  % in degs.

%% collect atmospheric data
    [mach,qbar] = ADC(vt,h);
    
elev = -(cm0+cma*alpha)/cmde;    % degrees

cl = cl0 + cla*alpha;
cd = dcdg + cd0 + cdcls*cl^2;

qs = qbar *s;       %% non-dimsionalizing factor

weight = mass * gd;
num_1 = qs*cl0-weight*cos(gam); num_2 = qs*cla;
den_1 = qs*(cd0+cdcls*cl0^2) + weight*sin(gam); 
den_2 = qs*2*cdcls*cl0*cla;
den_3 = qs*cdcls*cla^2;

[alpha,fval,flag] = fminsearch('feq',0);        %% alpha in degrees

cl = cl0 + cla*alpha;
cd = dcdg + cd0 + cdcls*cl^2;

lift = qs*cl; drag = qs*cd;
thrust = (drag+weight*sin(gam))/cos(alpha/rtod);

% res = thrust*sin(alpha/rtod)+lift-weight*cos(gam);        %% for debug
del_thr = thrust/(tstat+dtdv*vt);


elev = -((thrust*ze/(qs*cbar))+(cm0+dcmg+cma*alpha))/cmde;

% [del_thr elev alpha alpha/rtod]                   %% for debug
alp_r = alpha/rtod; theta_r = gam + alp_r;          %% alpha in radians

temp = [6,4,vt,alp_r,theta_r,0,h,0,del_thr,elev,.25,0];
name = input('name of output file? : ','s');
dlmwrite(name, temp);