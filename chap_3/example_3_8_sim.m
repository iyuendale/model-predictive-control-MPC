% constraints applied on the amplitude of control signals
num = 1;
den = [1 2*0.1*3 9];
[ac, bc, cc, dc] = tf2ss(num, den);
[Ap, Bp, Cp, Dp] = c2dm(ac, bc, cc, dc, 0.1);
%% augmented system
A = [Ap zeros(size(Ap, 1), size(Cp, 1));
	  Cp*Ap eye(size(Cp, 1))];
B = [Bp; Cp*Bp];
C = [zeros(size(Cp, 1), size(Ap, 1)) 1];
D = 0;
%% parameters
Q = C'*C; R = 0.3; Np = 46; x = [0.1 0.2 0.3]';
a = 0.7; N = 8; N_sim = 46;
%% DLQR
K = dlqr(A, B, Q, R);
y = []; kk = []; deltaU = []; u = 6; U = [];
for k = 10:56
	y = [y; k C*x];
	deltau = -K*x; u = deltau + u;
	x = A*x+B*deltau;
	deltaU = [deltaU; [k k+1]' [deltau deltau]'];
	U = [U; [k k+1]' [u u]'];
end
k = 10:56;
figure(1)
subplot 311, plot(y(:, 1), y(:, 2)), hold on
subplot 312, plot(deltaU(:, 1), deltaU(:, 2)), hold on
subplot 313, plot(U(:, 1),U(:, 2)), hold on 

%% DMPC - unconstrained
[Al, L0] = lagd(a, N);
L = zeros(N, N_sim);
for m = 1:Np
	L(:, m) = Al^(m-1)*L0;
end
% 
% % omega and psi
% omega = 0; psi = 0; R_L = R*eye(N);
% for m = 1:Np
% 	phi = 0;
% 	for j =  0:m-1
% 		phi = phi + (A^(m-1-j)*B*L(:, j+1)')';
% 	end
% 	omega = omega + phi*Q*phi';
% 	psi = psi + phi*Q*A^m;
% end
% omega = omega + R_L;

[omega, psi] = dmpc(A, B, a, N, Np, Q, R);
% optimal coefficients
x = [0.1 0.2 0.3]';
n = -(omega\psi)*x;
Kmpc = L(:, 1)'*(omega\psi);
y_dmpc = []; deltaU_dmpc = []; u = 6; U_dmpc = [];
for i = 10:56
	deltau = -Kmpc*x; u = deltau + u;
	y_dmpc = [y_dmpc; i C*x];
	x = A*x + B*deltau;
	deltaU_dmpc = [deltaU_dmpc; [i i+1]' [deltau deltau]'];
	U_dmpc = [U_dmpc; [i i+1]' [u u]'];
end
figure(1)
subplot 311, plot(y_dmpc(:, 1), y_dmpc(:, 2))
subplot 312, plot(deltaU_dmpc(:, 1), deltaU_dmpc(:, 2))
subplot 313, plot(U_dmpc(:, 1), U_dmpc(:, 2))

%% DMPC - constrained, constraints on both u and deltau
M_u = M_u_const(a, N, Np);
[M_du1, Lzerot] = Mdu(a, N, 1, Np);
M_du = [M_du1; -M_du1];
M = [M_u; M_du];
y_cons = []; x = [0.1 0.2 0.3]';
deltaU_cons = []; E = omega;
U_cons = [];
u_max = 4; u_min = 1.8; du_max = 0.25; du_min = -1;
u = 6;
for i = 10:56
	y_cons = [y_cons; i C*x];
	b =  [(u_max - u)*ones(Np, 1); (-u_min + u)*ones(Np, 1);
		       du_max*ones(Np, 1); -du_min*ones(Np, 1)];
	F = psi*x;
	deltau = QPhild(E, F, M, b); u = u + deltau(1);
	x = A*x + B*deltau(1);
	deltaU_cons = [deltaU_cons; [i i+1]' [deltau(1) deltau(1)]'];
	U_cons = [U_cons; [i i+1]' [u u]'];
end

figure(1)
subplot 311, plot(y_cons(:, 1), y_cons(:, 2),'k', 'linewidth', 1), grid on
legend 'DLQR' 'unconstrained DMPC' 'constrained DMPC'
title 'Output', ylabel 'y', xlabel 'k', axis([10 60 -0.5 1])
subplot 312, plot(deltaU_cons(:, 1), deltaU_cons(:, 2),'k', 'linewidth', 1)
yline(-1, 'k--'), yline(0.25, 'g--')
legend 'DLQR' 'unconstrained DMPC' 'constrained DMPC' ...
	'\Deltau_{min}' '\Deltau_{max}'
title 'Control increment', ylabel '\Deltau', xlabel 'k', grid on
subplot 313, plot(U_cons(:, 1), U_cons(:, 2),'k', 'linewidth', 1), grid on
yline(1.8, 'k--'), yline(4, 'g--')
legend 'DLQR' 'unconstrained DMPC' 'constrained DMPC' ...
	'u_{min}' 'u_{max}'
title 'Control variable', ylabel 'u', xlabel 'k'