% reduced weighted matrix R from 11 to 0.1
%% model of the system
am = 0.8; bm = 0.6; cm = 1;
% augmented system
A = [am 0; cm*am 1];
B = [bm; cm*bm]; C = [0 1];
%% given
ki = 10;    % initial time instance
Np = 16; Q = [0 0; 0 1]; R = 0.1; a = 0.6; N = 5;
%% Laguerre network parameters
x = [0.1 0.2]';
% Laguerre networks
[Al, L0] = lagd(a, N);
L = []; L(:, 1) = L0;
N_sim = 16;
for k = 2:N_sim
	L(:, k) = Al*L(:, k-1);
end

% optimal solution
R_L = R*eye(N); omega = 0; psi = 0;
for m = 1:Np
	phi = 0;
	for i = 0:m-1
		phi = phi + (A^(m-i-1)*B*L(:, i+1)')';
	end
	omega = omega + phi*Q*phi';
	psi = psi + phi*Q*A^m;
end
omega = omega + R_L;
% optimal coefficients
n = -inv(omega)*psi*x
buf = [10 C*x]; buf2 =[];
x1 = x;
%% DLQR
K = dlqr(A, B, Q, R);
disp(['DLQR gain: [', num2str(K), ']'])
buf3 = [10 C*x]; xlqr = x; buf4 = [];
sum = 0;        % sum of squared error
for k = 1:Np;
	% MPC
	deltau_mpc = n'*L(:, k);
	buf2 = [buf2; [k+ki-1 k+ki]' [deltau_mpc deltau_mpc]'];
	x1 = A*x1 + B*deltau_mpc;
	buf = [buf; k+ki C*x1];
	% DLQR
	deltau_lqr = -K*xlqr;
	xlqr = A*xlqr + B*deltau_lqr;
	buf3 = [buf3; k+ki C*xlqr];
	buf4 = [buf4; [k+ki-1 k+ki]' [deltau_lqr deltau_lqr]'];
	% sum of squared error
	sum = sum+ (deltau_lqr - deltau_mpc)^2;
end
figure(1)
subplot 211, yline(0, '--'), hold on, plot(buf(:, 1), buf(:, 2))
subplot 212, yline(0, '--'), hold on, plot(buf2(:, 1), buf2(:, 2))

subplot 211, plot(buf3(:, 1), buf3(:, 2)), hold off
subplot 212, plot(buf4(:, 1), buf4(:, 2)), hold off
subplot 211, legend('Set point', ['N = ', num2str(N)], 'DLQR'), title 'output'
subplot 212, legend('Set point', ['N = ', num2str(N)], 'DLQR'), title 'optimal control'

disp(['squared error for N = ', num2str(N), ' is: ', num2str(sum)])