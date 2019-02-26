function [MACH, QBAR] = ADC(Vt, h)

s_density = 2.377E-3; % Sea level density (slug/ft^3)
TFAC = 1.0 - 0.703e-5 * h;
T = 519.0*TFAC; % Temperature

if (h >= 35000.0)
    T = 390.0;
end
RHO = s_density * (TFAC^4.14); % Density
MACH = Vt/sqrt(1.4*1716.3*T); % Mach Number
QBAR = 0.5*RHO*Vt*Vt; % Dynamic Pressure
%C PS = 1715.0 * RHO * T 