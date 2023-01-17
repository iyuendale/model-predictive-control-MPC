% constraints applied on the difference of control signals
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
X = zeros(size(x, 1), N_sim);
X(:, 1) = x; kk = [];
U = zeros(size(B, 2), 50);
for k = 2:46
	u = -K*x;
	x = A*x+B*u;
	X(:, k) = x;
	U(:, k+8:k+9) = [u u];
	kk = [kk; [k+8:k+9]' [u u]'];
end
k = 10:55;
figure(1)
subplot 211, plot(k, X(3, :)), hold on
subplot 212, plot(kk(:, 1), kk(:, 2)), hold on

%% DMPC - unconstrained
[Al, L0] = lagd(a, N);
L = zeros(N, N_sim);
for m = 1:Np
	L(:, m) = Al^m*L0;
end
% omega and psi
omega = 0; psi = 0; R_L = R*eye(N);
for m = 1:Np
	phi = 0;
	for j =  0:m-1
		phi = phi + (A^(m-1-j)*B*L(:, j+1)')';
	end
	omega = omega + phi*Q*phi';
	psi = psi + phi*Q*A^m;
end
omega = omega + R_L;
% optimal coefficients
x = [0.1 0.2 0.3]';
n = -(omega\psi)*x;
Kmpc = L(:, 1)'*(omega\psi);
buf = []; buf2 = [];
for i = 10:60
	deltau = -Kmpc*x;
	buf = [buf; i C*x];
	x = A*x + B*deltau;
	buf2 = [buf2; [i i+1]' [deltau deltau]'];
end
figure(1)
subplot 211, plot(buf(:, 1), buf(:, 2))
title 'Output', ylabel 'y', xlabel 'k'
legend 'DLQR' 'unconstrained DMPC'
subplot 212, plot(buf2(:, 1), buf2(:, 2))
title 'Control increment', ylabel '\Deltau', xlabel 'k'
legend 'DLQR' 'unconstrained DMPC'
axis([10 60 -1.25 0.5])

%% DMPC - constrained
buf3 = []; x = [0.1 0.2 0.3]';
buf4 = [];
E = omega; Nc = 15;
%%% for constraints applied only on the ∆u(k)
% M_du = [-L(:, 1)'; L(:, 1)'];
% b_du = [1; 0.25];
%%% for costraints applied on values of Nc = 15 values of ∆u
[M_du1, Lzerot] = Mdu(a, N, 1, Nc);
b_du = [0.25*ones(Nc, 1); 1*ones(Nc, 1)];
M_du = [M_du1; -M_du1];  
% M_du = [L(:, 1:15)'; -L(:, 1:15)'];

for i = 10:60
	buf3 = [buf3; i C*x];
	F = psi*x;
	deltau = QPhild(E, F, M_du, b_du);
	x = A*x + B*deltau(1);
	buf4 = [buf4; [i i+1]' [deltau deltau]'];
end
figure(2)
subplot 211, plot(buf3(:, 1), buf3(:, 2))
title 'Output', ylabel 'y', xlabel 'k'
subplot 212, plot(buf4(1:end-1, 1), buf4(1:end-1, 2)), hold on
yline(0.25, 'k--'), yline(-1, 'k--')
title 'Control increment', ylabel '\Deltau', xlabel 'k'
axis([10 60 -1.25 0.5])