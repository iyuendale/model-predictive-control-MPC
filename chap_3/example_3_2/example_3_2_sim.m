clear all, clc
%% model of the system
am = 0.8; bm = 0.6; cm = 1;
% augmented system
A = [am 0; cm*am 1];
B = [bm; cm*bm]; C = [0 1];
%% given
ki = 10;    % initial time instance
x = [0.1 0.2]';
Np = 16; Q = [0 0; 0 1]; R = 1; a = 0.6; N = [1 2 3 4];
%% Laguerre network parameters
% for NN = 1:4
[Al, L0] = lagd(a, N(1))
L = [];
L(:, 1) = L0;
N_sim = 16;
for k = 2:N_sim
	L(:, k) = Al*L(:, k-1);
end

%% optimal solution
R_L = R*eye(N(1)); omega = 0; psi = 0;
for m = 1:Np
	phi = 0;
	for i = 0:m-1
		phi = phi + (A^(m-i-1)*B*L(:, i+1)')';
	end
	omega = omega + phi*Q*phi' + R_L;
	psi = psi + phi*Q*A^m;
end
% phi, omega, psi
n = -inv(omega)*psi*x
buf = []; buf2 =[];
for k = 1:Np;
	deltau = n'*L(:, k);
	x = A*x + B*deltau;
	buf = [buf; k C*x];
	buf2 = [buf2; [k k+1]' [deltau deltau]'];
end
subplot 211, plot(buf(:, 1), buf(:, 2), 'linewidth', N(1)), hold on
subplot 212, plot(buf2(:, 1), buf2(:, 2), 'linewidth', N(1)), hold on
% end
