%% 1
% state space model
tic
num = 0.01; den = [1 -0.6]
G = tf(num, den)
[Am, Bm, Cm, Dm] = tf2ss(num, den)
Bm = 0.01*Bm, Cm = 1
% augmented system
A = [Am zeros(size(Am, 1), size(Cm, 2)); Cm*Am eye(size(Cm,1))]
B = [Bm; Cm*Bm], C = [zeros(size(Cm, 1), size(Am,1)) Cm]
%% model predictive control system
Np = 16; Nc = 4; R = eye(Nc, Nc); r = 0.5;
Rs = r*ones(Np, 1);
F = []; phi = [];
for k = 1:Np
	F = [F; C*A^(k)];
	phi_col = [];
	for j = 1:Nc
		if k<j
			phi_col = [phi_col zeros(size(C,1), size(B, 2))];
		else
			phi_col = [phi_col C*A^(k-1)*B];
		end
	end
	phi = [phi; phi_col];
end
F, phi
Ky = inv(phi'*phi+R)*phi'*Rs;
Kmpc = inv(phi'*phi+R)*phi'*F;
%% 2
x = [0 0]'; t = 0; buf = [t x' 0 0]; u = 0; X= 0;
for k = 1:600
	deltaU = Ky-Kmpc*x;
	x = A*x + B*deltaU(1);   % augmented system
	u = u + deltaU(1);
	X = Am*X+Bm*u;   % plant
	buf = [buf; k x' deltaU(1) X]; 
end
t = buf(:, 1); x1 = buf(:, 2); x2 = buf(:, 3);
deltaU = buf(:, 4); X = buf(:, 5);
subplot 311, plot(t, x2)
subplot 312, plot(t, X)
subplot 313, plot(t, deltaU)
toc