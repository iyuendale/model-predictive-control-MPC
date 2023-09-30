% Predictive Control within One Optimization Window

% discrete system model
a = 0.8; b = 0.1; c = 1;

% augmented system
A = [a zeros(size(a, 1), size(c, 1)); c*a eye(size(c, 1))];
B = [b; c*b];
C = [zeros(size(c, 1), size(a, 1)) eye(size(c, 1))];

% components that form future output predictions
% Y = F*x(k) + phi*delta_U
Np = 10; Nc = 4;
F = zeros(Np*size(C, 1), size(A, 2));
phi = zeros(Np*size(C, 1), Nc*size(B, 2));
for i = 1:Np
	F(i, :) = C*A^i;
	for j = 1:Nc
		if i >= j
			phi(i, j) = C*A^(i-j)*B;
		end
	end
end

% controller parameters
phi_phi = phi'*phi
phi_F = phi'*F
RRs = ones(Np*size(C, 1), size(C, 1));
phi_RRs = phi'*RRs

%% optimal solution
r = 1; x = [0.1; 0.2]; rw1 = 0; u = 0; U = []; X = []; I = [];
xm = x(2);		% xm = y
R1 = rw1*eye(Nc);
figure(1)
for i = 10:14
	disp(['k = ', num2str(i)])
	delta_U1 = -(R1 + phi_phi)\(phi_F*x - phi_RRs*r)
	u = u + delta_U1(1)
	xm = x(2); X = [X; x'];
	xm_new = a*xm + b*u
	x = [xm_new-xm; xm_new]
	U = [U; u]; I = [I; i];
end
subplot(121), stairs(I, U, 'b')
title '(a) Optimal control', xlabel 'Sampling Instant'
ylabel 'Optimal Control', axis([9 15 1 8])
subplot(122), plot(I(1:end-1, 1), X(1:end-1, 1), 'b-*', I(1:end-1, 1), X(1:end-1, 2), 'k-o')
title '(b) State Variable', xlabel 'Sampling Instant'
ylabel 'State Variables', axis([9 15 0 1.25])
legend '\Deltax' 'y'
