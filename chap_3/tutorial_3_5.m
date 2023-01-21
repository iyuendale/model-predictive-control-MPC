% examples for MIMO system unconstrained control
%% the transfer function of a MIMO system
% there are 5 inputs and 4 outputs
g11 = tf(0.34, [0.85 1]); g12 = tf(0.21, [0.42 1]);
g13 = tf([0.25 0.5], [12 0.4 1]); g14 = 0; g15 = tf(6.46*[0.9 1], [0.07 0.3 1]);
g21 = tf(-0.41, [2.41 1]); g22 = tf(0.66, [1.51 1]);
g23 = tf(-0.3, [1.45 1]); g24 = 0; g25 = tf(-3.72, [0.8 1]);
g31 = tf(0.3, [2.54 1]); g32 = tf(0.49, [1.54 1]);
g33 = tf(-0.71, [1.35 1]); g34 = tf(-0.2, [2.72 1]);
g35 = tf(-4.71, [0.008 0.41 1]); g41 = 0;
g42 = 0; g43 = 0; g44 = 0;
g45 = tf(1.02, [0.07 0.31 1]);
G = [g12 g12 g13 g14 g15;
	  g22 g22 g23 g24 g25;
	  g32 g32 g33 g34 g35;
	  g42 g42 g43 g44 g45];
  
%% the minimal realization by using the 'ss' command
h = 1;
Gsmin = ss(G, 'min');
[Ac, Bc, Cc, Dc] = ssdata(Gsmin);
[Ap, Bp, Cp, Dp] = c2dm(Ac, Bc, Cc, Dc, h, 'zoh');
[m1, n1] = size(Cp); [n1, n_in] = size(Bp);
disp(['number of inputs: ', num2str(n_in)])
disp(['number of outputs: ', num2str(m1)])
disp(['number of states: ', num2str(n1)])

%% specify the the parametrs in Laguerre functions for each input
a = [0.5 0.5 0.5 0.5 0.5]; % a = 0.5 for each input
N = [10 10 10 10 10];  % N = 10 for each input
Np = 100;

%% augmented system of the plat with integrators and Q & R
% augmented system state equations
A = [Ap, zeros(n1, m1);
	  Cp*Ap eye(m1)];
 B = [Bp; Cp*Bp];
 C = [zeros(m1, n1) eye(m1)];
%  weighting matrices
Q = C'*C; R = 0.1*eye(n_in);
 %% generate omega and psi matrices by calling 'dmpc'
 [omega, psi] = dmpc(A, B, a, N, Np, Q, R);
 Lm = zeros(n_in, sum(N));
 [Al, L0] = lagd(a(1), N(1));
 L_m(1, 1:N(1)) = L0';   % for the first input
 In_s = 1;
 for jj= 2:n_in
	 [Al, L0] = lagd(a(jj), N(jj));
	 In_s = N(jj-1)+In_s;
	 In_e = In_s+N(jj)-1;
	 L_m(jj, In_s:In_e) = L0';
 end
 K = L_m*(omega\psi);
 A_closed = A - B*K;
%  eigenvalues of the closed loop
eig(A_closed);

%% state feedback using DLQR
[K_lqr, S, eig_dlqr] = dlqr(A, B, Q, R);
figure(1)
plot(eig_dlqr, 'ro'), hold on
plot(eig(A_closed), 'b*')
legend 'eigenvalues from DLQR' 'eigenvalues form MPC'
grid on

%% closed loop simulation
% initial conditions
y = zeros(m1, 1); u = zeros(n_in, 1); xm = zeros(n1, 1);
N_sim = 100;
% setpoints 
r1 = ones(1, N_sim+10);
r2 = zeros(1, N_sim+10);
r3 = zeros(1, N_sim+10);
r4 = zeros(1, N_sim+10);
sp = [r1; r2; r3; r4];
[M, Lzerot] = Mdu(a, N, n_in, 1);
[u1, y1, deltau1, k] = simuuc(xm, u, y, sp, Ap, Bp, Cp, ...
								N_sim, omega, psi, Lzerot);
subplot 211, plot(k, y1)
subplot 212, plot(k, u1)