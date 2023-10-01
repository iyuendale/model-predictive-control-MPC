% How to implement a predictive control system

Am = [1 1; 0 1]; Bm = [0.5; 1];
Cm = [1 0]; Dm = 0;
Np = 20; Nc = 4; rs = 1;

% Controller gains
[phi_phi,phi_F,phi_Rs,AA,BB,CC] = mpcgains(Am,Bm,Cm,Np,Nc,rs);
% Initial values
xm = [0; 0]; u = 0; y = 0;
x = zeros(3, 1);
rw = 0.1;   % change the value to 2 to see the response, slower
N_sim = 100;
R = rw*eye(Nc);

Ky = [1 0 0 0]*inv(phi_phi + R)*phi_Rs;
Kmpc = [1 0 0 0]*inv(phi_phi + R)*phi_F;
U = [0 u]; Y = [0 y];

for i = 1:N_sim
	deltau = Ky - Kmpc*x;
	u = u + deltau;
	xm_new = Am*xm + Bm*u;
	y = Cm*xm_new;
	x = [(xm_new-xm); y]; xm = xm_new;
	U = [U; i u]; Y = [Y; i y];
end

subplot 211
plot(Y(:, 1), Y(:, 2)), legend 'Output'
xlabel 'Sample instant', axis([0 100 0 1.25])
subplot 212
stairs(U(:, 1), U(:, 2)), legend 'Control input'
xlabel 'Sample instant', axis([0 100 -0.5 0.8])