% closed-loop feedback gain matrices

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
phi_phi = phi'*phi;
phi_F = phi'*F;
RRs = ones(Np*size(C, 1), size(C, 1));
phi_RRs = phi'*RRs;

%% closed-loop feedback gain matrices
r = 1; rw1 = 0; rw2 = 10;
disp(['r_w = ', num2str(rw1)])
R1 = rw1*eye(Nc)
Ky = [1 0 0 0]*inv(phi_phi+R1)*phi_RRs
Kmpc = [1 0 0 0]*inv(phi_phi+R1)*phi_F
AA = A - B*Kmpc
eigenvalue = eig(AA)
disp(['r_w = ', num2str(rw2)])
R2 = rw2*eye(Nc);
Ky = [1 0 0 0]*inv(phi_phi+R2)*phi_RRs
Kmpc = [1 0 0 0]*inv(phi_phi+R2)*phi_F
AA = A - B*Kmpc
eigenvalue = eig(AA)