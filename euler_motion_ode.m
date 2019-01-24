function dwdt = euler_motion_ode(t, w, J)
dwdt = zeros(3,1);
dwdt(1) = ((J(2) - J(3)) * w(2) * w(3)) / J(1);
dwdt(2) = ((J(3) - J(1)) * w(3) * w(1)) / J(2);
dwdt(3) = ((J(1) - J(2)) * w(1) * w(2)) / J(3);
%Omega = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
%dwdt = -inv(J)*(Omega*J*w);
