% producing stable closed-loop predictive control system from
% unstable system with one pole and zero are outside of the 
% unit circle
%% discrete time model of the system
numd = [0.3 -0.33]; dend = conv([1 -0.3], [1 -1.5]);
[am, bm, cm, dm] = tf2ss(numd, dend);
% augmented system
A = [am zeros(size(am, 1), size(cm, 1)); cm*am eye(size(cm, 1))];
B = [bm; cm*bm];
C = [zeros(size(cm, 1), size(am, 1)) eye(size(cm, 1))]; D = 0;

%% Control parameters
% Laguerre parameters
a = 1/1.5;		% reciprocal of the unstable pole
N = 3;
[Al, L0] = lagd(a, N);
% predictive control parameters
% let Np (prediction horizon) ranges from 18 to 40 with step of 2
NP = 18:2:40;

figure(1)
subplot 311, title output, hold on
subplot 312, title 'control variable', hold on
subplot 313, title 'control increment', hold on
for iii = 1:size(NP, 2)
	Np = NP(iii);			% change the value of to 18 and 38
	Q = C'*C; R = 0.1;
	% omega and psi
	[omega, psi] = dmpc(A, B, a, N, Np, Q, R);
	omega
	Kmpc = L0'*(omega\psi);
	% initials
	x = zeros(size(A, 1), 1); buf = []; u = 0;
	xm = zeros(size(am, 1), 1); xm_old = xm;
	for k = 1:100
		deltau = -Kmpc*x;
		u = u+deltau;
		xm = am*xm+bm*u; y = cm*xm;
		x = [xm-xm_old; y-1]; xm_old = xm;
		buf = [buf; k deltau u y];
	end
	figure(1)
	subplot 311, plot(buf(:, 1), buf(:, 4)), hold on
	subplot 312, plot(buf(:, 1), buf(:, 3)), hold on
	subplot 313, plot(buf(:, 1), buf(:, 2)), hold on
% 	pause(1)
end
% figure(1)
% subplot 311, title output
% subplot 312, title 'control variable'
% subplot 313, title 'control increment'

% As Np increases, the system reaches the set point with faster
% speed and is more stable.
% But the values in Hessian matrix, omega, become larger and 
% only the first element will be there with large value for the
% largest values of Np.