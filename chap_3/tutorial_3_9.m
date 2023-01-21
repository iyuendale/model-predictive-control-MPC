% clear all, clc
% continous system
T = 1; K = 1;
num = K; den = [T 1 0];
% discretize and state space model
[numd, dend] = c2dm(num, den, 0.1);
ad = [-dend(2:3) numd(3); 1 0 0; 0 0 0];
bd = [numd(2) 0 1]'; cd = [1 0 0];
% augmented system
A = [ad zeros(size(ad, 1), size(cd, 1));
	  cd*ad eye(size(cd, 1))];
B = [bd; cd*bd];
C = [zeros(size(cd, 1), size(ad, 1)) eye(size(cd, 1))];
% design parameters
Q = C'*C; R = 0.1; Np = 46;
% Laguerre parameters
N = 1; a = 0.5;
[omega, psi] = dmpc(A, B, a, N, Np, Q, R);
[Al, L0] = lagd(a, N);
Kmpc = L0'*(omega\psi);
A_closed = A-B*Kmpc;
% closed loop poles
eig_closed = eig(A_closed);
% parameters for control simulation
N_sim = 170;
u_min= -0.3; u_max= 0.2; deltau_min= -0.1; deltau_max= 0.1;
u = 0; r = 1; xm = [0 0 0]'; xm_old = xm; y = 0;
xf = [xm-xm_old; y-r];
buf = []; buf2 = [];

%% unconstrained control
for k = 1:N_sim
	deltau = -Kmpc*xf;
	if (deltau > deltau_max) deltau = deltau_max; end
	if (deltau < deltau_min) deltau = deltau_min; end
	u = u + deltau;
	if (u > u_max) u = u_max; end
	if (u < u_min) u = u_min; end
	xm = ad*xm + bd*u; y = cd*xm;
	xf = [xm-xm_old; y-r]; xm_old = xm;
	buf = [buf; k y];
	buf2 = [buf2; [k k+1]' [deltau deltau]' [u u]'];
end
subplot 211, plot(buf2(:, 1), buf2(:, 3)), hold on
yline(u_min), yline(u_max)
subplot 212, plot(buf(:, 1), buf(:, 2))