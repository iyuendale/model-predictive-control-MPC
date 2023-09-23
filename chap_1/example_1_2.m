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
r = 1; x = [0.1; 0.2]; rw1 = 0; rw2 = 10;
R1 = rw1*eye(Nc); R2 = rw2*eye(Nc);
delta_U1 = -(R1 + phi_phi)\(phi_F*x - phi_RRs*r)
delta_U2 = -(R2 + phi_phi)\(phi_F*x - phi_RRs*r)

%% comparison of optimal solutions
Y1 = F*x + phi*delta_U1;
Y2 = F*x + phi*delta_U2;
y1 = [C*x; Y1]; y2 = [C*x; Y2];
delta_xm1 = [x(1)];
delta_xm2 = [x(1)];
for i = 2:11
	delta_xm1 = [delta_xm1; y1(i) - y1(i-1)];
	delta_xm2 = [delta_xm2; y2(i) - y2(i-1)];
end
t = 10:1:20;
subplot 211
plot(t, y1, '-o', t, delta_xm1, '-*'), axis([10 20 -.2 1.2])
title 'State variables with no weight on \Deltau'
legend 'y' '\Deltax_m'
xlabel 'Sampling instant', ylabel 'Response'

subplot 212
plot(t, y2, '-o', t, delta_xm2, '-*'), axis([10 20 0 0.8])
title 'State variables with weight on \Deltau'
legend 'y' '\Deltax_m'
xlabel 'Sampling instant', ylabel 'Response'
