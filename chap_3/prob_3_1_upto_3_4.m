%% problem 3.1
% Discretization
num = 1; den = conv([10 1], conv([10 1], [10 1]));
sys = tf(num, den, 'Variable', 's', 'IODelay', 3)
[num1, den1] = tfdata(sys);
[numd, dend] = c2dm(num1, den1, 1);
sysd = tf(numd, dend, 1)
am = [-dend(2) -dend(3) -dend(4) numd(3) numd(4);
	     1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0; 0 0 0 1 0];
bm = [numd(2) 0 0 1 0]';
cm = [1 0 0 0 0]; dm = 0;
AA = [am zeros(size(am, 1), size(cm, 1)); cm*am eye(size(cm, 1))];
BB = [bm; cm*bm]; CC = [zeros(1, size(am, 1)) 1]; DD = 0;
% Impulse response
H = dimpulse(am, bm, cm, dm);

% Laguerre parameters
a = 0.9; N = 3;
[Al, L0] = lagd(a, N);
N_sim = size(H, 1);
L= [];
for k = 1:N_sim
	L(:, k) = Al^(k-1)*L0;
end

% coefficients
c = L(:, :)*H;

% impulse response of the approximated model
H_model = c'*L;
figure(1), plot([H, H_model'], 'linewidth', 2)
legend 'Plant' 'Laguerre Model', xlabel 'sampling instant', ylabel 'response'
title 'Impulse response approximation'
%% problem 3.2
z = tf('z', 1);
laguerre_fun = 1;
for i = 1:N-1
	laguerre_fun = [laguerre_fun; ((z^-1-a)/(1-a*z^-1))^i];
end

GA = (sqrt(1-a^2)/(1-a*z^-1))*c'*laguerre_fun;
[numd2, dend2] = tfdata(GA); kk = 1:1:size(H_model, 2);
figure(2), dimpulse(numd2, dend2, kk); hold on
xm = L; Cm = c'; Am = Al; Bm = L0;
ym = Cm*xm;
plot(kk, ym, 'linewidth', 2)
legend 'G_A(z)' 'Laguerre model'
title 'State-Space Relization'
% augmented system
A = [Am zeros(size(Am, 1), size(Cm, 1)); Cm*Am eye(size(Cm, 1))];
B = [Bm; Cm*Bm]; C = [zeros(1, size(Am, 1)) 1]; D = 0;

%% problem 3.3
Np = 30; N = 3; a = 0.5; 
Q = C'*C; R = 1; r = 1; d = 1;
[Al, L0] = lagd(a, N);
L = [];
for k = 1:Np
	L = [L Al^(k-1)*L0];
end
[omega, psi] = dmpc(A, B, a, N, Np, Q, R);
Kmpc = -L0'*(omega\psi);
% using the Laguerre model G_A(z)
x = zeros(size(A, 1), 1); buf = []; xm = zeros(size(Am, 1), 1);
xm_old = xm; u = 0;
for k = 1:29
	deltau = Kmpc*x;
	u = u+deltau;
	xm = Am*xm + Bm*u;
	x = [xm-xm_old; Cm*xm-r]; xm_old = xm;
	buf = [buf; deltau u Cm*xm];
end

for k = 30:100
	deltau = Kmpc*x;
	u = u+deltau;
	xm = Am*xm + Bm*(u-d);
	x = [xm-xm_old; Cm*xm-r]; xm_old = xm;
	buf = [buf; deltau u Cm*xm];
end
figure(3)
subplot 311, plot(buf(:, 3)), hold on
subplot 312, stairs(buf(:, 1)), hold on
subplot 313, stairs(buf(:, 2)), hold on

% using the discret time model G(z)
Q2 = CC'*CC;
[omega2, psi2] = dmpc(AA, BB, a, N, Np, Q2, R);
Kmpc2 = -L0'*(omega2\psi2);
x = zeros(size(am, 1)+1, 1); buf = []; xm = zeros(size(am, 1), 1);
xm_old = xm; u = 0;
for k = 1:29
	deltau = Kmpc2*x;
	u = u+deltau;
	xm = am*xm + bm*u;
	x = [xm-xm_old; cm*xm-r]; xm_old = xm;
	buf = [buf; deltau u cm*xm];
end

for k = 30:100
	deltau = Kmpc2*x;
	u = u+deltau;
	xm = am*xm + bm*(u-d);
	x = [xm-xm_old; cm*xm-r]; xm_old = xm;
	buf = [buf; deltau u cm*xm];
end
figure(3)
subplot 311, plot(buf(:, 3))
legend 'Laguerre model' 'Discret time model'
title 'Output/Response', xlabel 'sample instant-k', ylabel 'y[k]'
subplot 312, stairs(buf(:, 1))
legend 'Laguerre model' 'Discret time model'
title 'Control Increment', xlabel 'sample instant-k', ylabel '\Deltau[k]'
subplot 313, stairs(buf(:, 2))
legend 'Laguerre model' 'Discret time model'
title 'Control Signal', xlabel 'sample instant-k', ylabel 'u[k]'

%% problem 3.4
% change the value of N to 8 in the above problems to see the 
% performance changes