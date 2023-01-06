%% the discrete time model dynamics
am = [0.9048 0; 0.0952 1]; bm = [0.0952; 0.0048];
cm = [0 1]; dm = 0;

% augmented system
A = [am zeros(size(am, 1), 1); cm*am eye(size(cm, 1))]
B = [bm; cm*bm]
C = [zeros(size(cm, 1), size(am, 1)) eye(size(cm, 1))], D = 0

%% given parameters
Np = 20; Nc = 4; R = eye(Nc); rs = 1; Rs = rs*ones(Np, 1);
% observer poles
Poles = [0.1 0.2 0.3];
% observer gain
Kob = place(A', C', Poles)'

%% MPC controller
[F, phi]= F_phi(A, B, C, Np, Nc);
% constraints
M = [-1 zeros(1, Nc-1); 1 zeros(1, Nc-1); -1 zeros(1, Nc-1); 1 zeros(1, Nc-1)];
% initials
x = [0 0 0]'; x_est = [0 0 0]'; u = 0; buf = []; buf2 = [];
E = phi'*phi+R;

for k = 1:100
	FF = phi'*(F*x_est - Rs);
	b = [u; 0.6-u; 0.2; 0.2];
	deltaU = QPhild(E, FF, M, b); u = u + deltaU(1);
	x_est = A*x_est + B*deltaU(1) + Kob*(C*x - C*x_est);
	x = A*x +B*deltaU(1);
	buf = [buf; k, C*x_est C*x deltaU(1)];
	buf2 = [buf2; [k k+1]' [u u]'];
end

k = buf(:, 1); y_est = buf(:, 2); y = buf(:, 3); deltau = buf(:, 4);
k2 = buf2(1:end-1, 1); u = buf2(1:end-1, 2);
subplot 211, plot(k, y)
subplot 212, plot(k2, u)