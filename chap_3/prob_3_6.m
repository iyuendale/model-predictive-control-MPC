% design predictive control system for a vessel with two tanks 
% connected in series
% dynamic model of the system
alpha1 = 2; alpha2 = 1; beta = 1;
acm = [-beta/alpha1 beta/alpha1; beta/alpha2 -beta/alpha2];
bcm = [1/alpha1 0; 0 1/alpha2]; ccm = [1 0; 0 1]; dcm = zeros(2);
eigs_of_sys = eig(acm) % it has pole at 0, so integrator is inside it
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
a_values = [0 0.5 0.9]; N = [3 3];
for ai = 1:size(a_values, 2)
	a = [a_values(ai) a_values(ai)];
	L0_mimo = zeros(sum(N), size(a, 2));
	j = 0;
	for i = 1:size(a, 2)
		[Al, L0] = lagd(a(i), N(i));
		L0_mimo(j+1:N(i)+j, i) = L0;
		j = N(i);
	end
	
	% predictive control parameters
	Q = C'*C; R = 0.1*eye(size(C, 1)); Np = 20;
	[omega, psi] = dmpc(A, B, a, N, Np, Q, R);
	Kmpc = L0_mimo'*(omega\psi);
	
	x = zeros(size(A, 1), 1); xm = zeros(size(am, 1), 1); xm_old = xm;
	buf = []; u = 0;
	for k = 1:30
		deltau = -Kmpc*x;
		u = u + deltau;
		xm = am*xm + bm*u; y = cm*xm;
		x = [xm-xm_old; y-[1 1]']; xm_old = xm;
		buf = [buf; k deltau' u' y'];
	end
	
	figure(1)
	subplot 321, plot(buf(:, 1), buf(:, 6)), hold on
	subplot 322, plot(buf(:, 1), buf(:, 7)), hold on
	subplot 323, stairs(buf(:, 1), buf(:, 2)), hold on
	subplot 324, stairs(buf(:, 1), buf(:, 3)), hold on
	subplot 325, stairs(buf(:, 1), buf(:, 4)), hold on
	subplot 326, stairs(buf(:, 1), buf(:, 5)), hold on
end

figure(1)
subplot 321, legend 'a= 0' 'a= 0.5' 'a= 0.9'
title 'y_1'
subplot 322, legend 'a= 0' 'a= 0.5' 'a= 0.9'
title 'y_2'
subplot 323, legend 'a= 0' 'a= 0.5' 'a= 0.9'
title '\Deltau_1'
subplot 324, legend 'a= 0' 'a= 0.5' 'a= 0.9'
title '\Deltau_2'
subplot 325, legend 'a= 0' 'a= 0.5' 'a= 0.9'
title 'u_1'
subplot 326, legend 'a= 0' 'a= 0.5' 'a= 0.9'
title 'u_2'

% If the need is to have no overshoot in the response, use Laguerre
% scaling parameter of a = 0.9