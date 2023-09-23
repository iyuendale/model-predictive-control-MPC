%% the transfer function of a MIMO system
% there are 5 inputs and 4 outputs
g11 = tf([0.21048 0.00245], [1 0.302902 0.066775 0.002186]);
g12 = tf([-0.001313 0.000548 -0.000052], [1 0.210391 0.105228 0.00777 0.000854]);
g21 = tf([0.000976 -0.000226], [1 0.422036 0.091833 0.003434]);
g22 = tf(-0.00017, [1 0.060324 0.006836]);
G = [g11 g12; g21 g22];
  
%% 1
dt = 3;
Gsmin = ss(G, 'min');
[ac, bc, cc, dc] = ssdata(Gsmin);
[Ap, Bp, Cp, Dp] = c2dm(ac, bc, cc, dc, dt, 'zoh');
[n_out, ~] = size(Cp); [n_state, n_in] = size(Bp);
disp(['number of inputs: ', num2str(n_in)])
disp(['number of outputs: ', num2str(n_out)])
disp(['number of states: ', num2str(n_state)])

%% 2
% specify the the parametrs in Laguerre functions for each input
a = [0.6 0.6]; % a = 0.6 for each input
N = [4 4];  % N = 4for each input
Np = 28;

% augmented system of the plat with integrators and Q & R
% augmented system state equations
A = [Ap, zeros(n_state, n_out);
	  Cp*Ap eye(n_out)];
B = [Bp; Cp*Bp];
C = [zeros(n_out, n_state) eye(n_out)];
%  weighting matrices
Q = C'*C; R = 6*eye(n_in);
% generate omega and psi matrices
[omega, psi] = dmpc(A, B, a, N, Np, Q, R);
Lm_0 = zeros(sum(N), n_in);
j = 0;
for i = 1:size(a, 2)
	[Al, L0] = lagd(a(i), N(i));
	Lm_0(j+1:j+N(i), i) = L0;
	j = j+N(i);
end
Kmpc = Lm_0'*(omega\psi)   % closed loop control gain
A_closed = A - B*Kmpc;
eigen_closed = eig(A_closed);
figure(1), plot(eigen_closed, 'b*')          % eigenvalues of closed loop

% % control the system
% xp = zeros(size(Ap, 1), 1); xp_old = xp; buf = [];
% x = [xp-xp_old; Cp*xp]; u = [0; 0]; r = [1; 1]; d = r;
% 
% for k = 1:100
% 	deltau = -Kmpc*x;
% 	u = u+deltau;
% 	xp = Ap*xp + Bp*(u- d); y = Cp*xp;
% 	x = [xp-xp_old; y-r]; xp_old = xp;
% 	buf = [buf; k deltau' u' y'];
% end
% figure(2)
% subplot 321, plot(buf(:, 6)), hold on, title 'y_1'
% subplot 322, plot(buf(:, 7)), hold on, title 'y_2'
% subplot 323, stairs(buf(:, 4)), hold on, title 'u_1'
% subplot 324, stairs(buf(:, 5)), hold on, title 'u_2'
% subplot 325, stairs(buf(:, 2)), hold on, title '\Deltau_1'
% subplot 326, stairs(buf(:, 3)), hold on, title '\Deltau_2'

%% 3
Qob = eye(n_state); Rob = 0.1*eye(n_out);
Kob = dlqr(Ap', Cp', Qob, Rob)'
%% 4
xp = zeros(size(Ap, 1), 1); xp_old = xp; buf = [];
x = [xp-xp_old; Cp*xp]; u = [0; 0]; r = [1; 1]; d = r;
xp_est = zeros(size(Ap, 1), 1); xp_old_est = xp_est;
x_est = [xp_est-xp_old_est; Cp*xp_est]; r = [3 0]';
uss = [300 2000]'; yss = [639 42]';
for k = 1:100
	deltau = -Kmpc*x_est;
	u = u+deltau;
	xp_est = Ap*xp_est + Bp*u + Kob*Cp*(xp - xp_est);
	x_est = [xp_est-xp_old_est; Cp*xp_est-r];  xp_old_est = xp_est;
	xp = Ap*xp + Bp*u; y = Cp*xp;
	x = [xp-xp_old; y-r]; xp_old = xp;
	buf = [buf; k deltau' (u'+uss') (y'+yss')];
end
figure(3)
subplot 321, plot(buf(:, 6)), hold on, title 'y_1'
subplot 322, plot(buf(:, 7)), hold on, title 'y_2'
subplot 323, stairs(buf(:, 4)), hold on, title 'u_1'
subplot 324, stairs(buf(:, 5)), hold on, title 'u_2'
subplot 325, stairs(buf(:, 2)), hold on, title '\Deltau_1'
subplot 326, stairs(buf(:, 3)), hold on, title '\Deltau_2'
