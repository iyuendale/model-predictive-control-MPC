%% continuous model to discrete model conversion
num = [1 -3];
den = [1 2*0.3*3 3^2];
[a, b, c, d] = tf2ss(num, den);
% discretize
[am, bm, cm, dm] = c2dm(a, b, c, d, 0.1);
% augmented system
A = [am zeros(size(am, 1), 1); cm*am eye(size(cm, 1))];
B = [bm; cm*bm]; C = [0 0 1]; D = 0;
%% Given
R = 0.3; Q = C'*C; Np = 36; ki = 10; x = [0.1 0.2 0.3]';
% DLQR control gain and closed loop poles
Klqr = dlqr(A, B, Q, R)
eigenvalues_lqr = eig(A - B*Klqr)
% Laguerre parameters
N = 8; a_values = [0 0.4 0.6 0.8]';
% loop through the values of 'a'
for j = 1:size(a_values)
	a = a_values(j);
	[Al, L0] = lagd(a, N);
	L = []; L(:, 1) = L0;
	N_sim = 36;
	for k = 2:N_sim
		L(:, k) = Al*L(:, k-1);
	end
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
	n = -inv(omega)*psi*x;
	% optimal control gain
	disp(['feedback gain and closed loop poles for a = ', num2str(a)])
	Kmpc = L(:, 1)'*inv(omega)*psi
	eigenvalues = eig(A-B*Kmpc)
end