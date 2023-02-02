% design predictive control system for a vessel with two tanks 
% connected in series
% dynamic model of the system
clear all
alpha1 = 2; alpha2 = 1; beta = 1;
acm = [-beta/alpha1 beta/alpha1; beta/alpha2 -beta/alpha2];
bcm = [1/alpha1 0; 0 1/alpha2]; ccm = [1 0; 0 1]; dcm = zeros(2);

% discretize with delta_t = 0.1
[am, bm, cm, dm] = c2dm(acm, bcm, ccm, dcm, 0.1);
% augmented system
A = [am zeros(size(am, 1), size(cm, 1));
	  cm*am eye(size(cm, 1))];
B = [bm; cm*bm];
C = [zeros(size(cm, 1), size(am, 1)) eye(size(cm, 1))];
D = dm;
% Laguerre parameters
% Since the laguerre parameters are equal: a1 = a2 & N1 = N2
a = [0.8 0.8]; N = [3 3];
%% unconstrained
L0_mimo = zeros(sum(N), size(a, 2)); j = 0;
for i = 1:size(a, 2)
	[Al, L0] = lagd(a(i), N(i));
	L0_mimo(j+1:N(i)+j, i) = L0;
	j = N(i);
end

% predictive control parameters
Q = C'*C; R = 0.1*eye(size(C, 1)); Np = 20; N_sim = 80;
[omega, psi] = dmpc(A, B, a, N, Np, Q, R);
Kmpc = L0_mimo'*(omega\psi)
% Set points
sp = [ones(1, N_sim/2) 2*ones(1, N_sim/2); 0.5*ones(1, N_sim)];
x = zeros(size(A, 1), 1); xm = zeros(size(am, 1), 1); xm_old = xm;
buf = []; u = 0;
for k = 1:N_sim
	deltau = -Kmpc*x;
	u = u + deltau;
	xm = am*xm + bm*u; y = cm*xm;
	x = [xm-xm_old; y-sp(:, k)]; xm_old = xm;
	buf = [buf; k deltau' u' y'];
end

figure(1)
subplot 321, plot(buf(:, 1), buf(:, 6)), hold on
subplot 322, plot(buf(:, 1), buf(:, 7)), hold on
subplot 323, stairs(buf(:, 1), buf(:, 2)), hold on
subplot 324, stairs(buf(:, 1), buf(:, 3)), hold on
subplot 325, stairs(buf(:, 1), buf(:, 4)), hold on
subplot 326, stairs(buf(:, 1), buf(:, 5)), hold on

%% constrained
Nc = N_sim;    % number of constraints that will be imposed
deltau1_min = -0.3; deltau1_max = 0.3;
deltau2_min = -0.1; deltau2_max = 0.1;
u1_min = 0; u1_max = 2;
u2_min = -2; u2_max = 0.5;
[Al_1, L0_1] = lagd(a(1), N(1));
[Al_2, L0_2] = lagd(a(2), N(2));
L_1= []; L_2= [];
for i = 1:Nc
	L_1 = [L_1; (Al_1^(i-1)*L0_1)'];
	L_2 = [L_2; (Al_2^(i-1)*L0_2)'];
end
E = omega;
% constraint matrices on u and deltau
M_du = zeros(4*Nc, sum(N));
M_du(1:2*Nc, 1:N(1)) = [-L_1(1:Nc, :); L_1(1:Nc, :)];
M_du(2*Nc+1:4*Nc, N(1)+1:sum(N)) = [-L_2(1:Nc, :); L_2(1:Nc, :)];

M_u = zeros(4*Nc, sum(N));
M_u(1, 1:N(1)) = -L_1(1, :); 
M_u(1+Nc, 1:N(1)) = L_1(1, :);
M_u(1+2*Nc, 1+N(1):end) = -L_2(1, :);
M_u(1+3*Nc, 1+N(1):end) = L_2(1, :);
	
for k = 2:Nc
	M_u(k, 1:N(1)) = -sum(L_1(1:k, :));
	M_u(k+Nc, 1:N(1)) = sum(L_1(1:k, :));
	M_u(k+2*Nc, 1+N(1):end) = -sum(L_2(1:k, :));
	M_u(k+3*Nc, 1+N(1):end) = sum(L_2(1:k, :));
end
M = [M_du; M_u];

gamma_du = [-deltau1_min*ones(Nc, 1); deltau1_max*ones(Nc, 1);
	                -deltau2_min*ones(Nc, 1); deltau2_max*ones(Nc, 1)];

x = zeros(size(A, 1), 1); xm = zeros(size(am, 1), 1); xm_old = xm;
u = [0; 0]; buf = [];
for k = 1:N_sim
	F = psi*x;
	gamma_u = [(-u1_min+u(1))*ones(Nc, 1);
				   (u1_max-u(1))*ones(Nc, 1);
	                     (-u2_min+u(2))*ones(Nc, 1);
				   (u2_max-u(2))*ones(Nc, 1)];
	gamma = [gamma_du; gamma_u];
	eta = QPhild(E, F, M, gamma);
 	deltau = [L0_1' zeros(1, N(2)); zeros(1, N(1)) L0_2']*eta;
	u = u+deltau;
	xm = am*xm+bm*u; y = cm*xm;
	x = [xm-xm_old; y-sp(:, k)]; xm_old = xm;
	buf = [buf; k deltau' u' y'];
end
figure(1)
subplot 321, plot(buf(:, 1), buf(:, 6)), hold on
subplot 322, plot(buf(:, 1), buf(:, 7)), hold on
subplot 323, stairs(buf(:, 1), buf(:, 2)), hold on
subplot 324, stairs(buf(:, 1), buf(:, 3)), hold on
subplot 325, stairs(buf(:, 1), buf(:, 4)), hold on
subplot 326, stairs(buf(:, 1), buf(:, 5)), hold on

figure(1)
subplot 321, legend 'unconstrained' 'constrained'
title 'y_1'
subplot 322, legend 'unconstrained' 'constrained'
title 'y_2'
subplot 323, yline(deltau1_min, 'k--'), yline(deltau1_max, 'k--')
title '\Deltau_1', legend 'unconstrained' 'constrained'
subplot 324, yline(deltau2_min, 'k--'), yline(deltau2_max, 'k--')
title '\Deltau_2', legend 'unconstrained' 'constrained'
subplot 325, yline(u1_min, 'k--'), yline(u1_max, 'k--'), title 'u_1'
legend 'unconstrained' 'constrained', axis([0 N_sim -0.5 4.5])
subplot 326, yline(u2_min, 'k--'), yline(u2_max, 'k--'), title 'u_2'
legend 'unconstrained' 'constrained', axis([0 N_sim -2.5 1.5])
