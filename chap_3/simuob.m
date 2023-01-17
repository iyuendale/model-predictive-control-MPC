%% tutorial 3.6
% produces the the closed-loop simulation for control systems
% with observer in the loop
% feedback is based on the estimation of x_hat
function [u1, y1, deltau1, k] = simuob(xm, u, y, sp, Ap, Bp, ...
				Cp, A, B, C, N_sim, omega, psi, K_ob, Lzerot)
% closed loop simulation without constraints
[ny, n] = size(C); [n, nu] = size(B);
x_hat = zeros(n, 1);
% u1 = zeros(nu, N_sim);
for kk = 1:N_sim
	Xsp = [zeros(n-ny, 1); sp(:, kk)];   % setpoints signals on y
	eta = -(omega\psi)*(x_hat - Xsp);
	deltau = Lzerot*eta;
	u = u+deltau;     % update u
	deltau1(:, kk) = deltau;
	u1(1:nu, kk) = u;       % keep u
	y1(1:ny, kk) = y;       % keep y
	x_hat = A*x_hat + K_ob*(y - C*x_hat) + B*deltau;
	                                        % u and y to generate x_hat(k+1)
	%%% plant simulation
	xm = Ap*xm + Bp*u;       % calculates xm(k+1)
	y = Cp*xm;                      % calculates y(k+1)
end
k = 0:(N_sim-1);