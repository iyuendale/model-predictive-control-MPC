%% tutorial 3.4
% Predictive unconstrained control system simulation
% observer is not used
% this function produces multivariable closed loop simulation
% for unconstrained control systems
% number of rows of set points is equal to number of otputs and
% number of columns is greater than simulation time(N_sim)
function [u1, y1, deltau1, k]= simuuc(xm, u, y, sp, Ap, Bp, ...
								Cp, N_sim, omega, psi, Lzerot)
[m1, n1] = size(Cp); [n1, n_in] = size(Bp);
Xf = [xm; (y-sp(:, 1))];        % initialize feedback state variable
for kk = 1:N_sim
	eta = -(omega\psi)*Xf;
	deltau = Lzerot*eta;
	u = u+deltau;
	deltau1(:, kk) = deltau;
	u1(1:n_in, kk) = u;
	y1(1:m1, kk) = y;
	% plant simulation
	xm_old = xm;
	xm = Ap*xm + Bp*u;      % calculates xm(k+1)
	y = Cp*xm;                     % calculates y(k+1)
	% updating feedback state variable
	Xf = [xm-xm_old; y-sp(:, kk+1)];
end
k = 0:N_sim-1;