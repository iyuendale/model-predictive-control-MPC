%% Examine sensitivity issues in the selection of design parameters

omega = 10;
num = omega^2;
den = [1 0.1*omega omega^2];
% transfer function to continous SS representation
[Ac, Bc, Cc, Dc] = tf2ss(num, den)
% discretization
[Am, Bm, Cm, Dm] = c2dm(Ac, Bc, Cc, Dc, 0.01)
% designing augmented system
A = [Am zeros(size(Am, 1), size(Cm, 1));
	  Cm*Am eye(size(Cm, 1))]
B = [Bm; Cm*Bm]
C = [zeros(size(Cm, 1), size(Am, 1)) eye(size(Cm, 1))]

Nc = 3; rw = 0.5; x = [0.1 0.2 0.3]';
R = rw*eye(Nc); rs = 1;

%% if Np = 20
Np = 20;
disp(['Np = ', num2str(Np)]);
[F,phi,RRs] = IEmpcgains2(A,B,C,Np,Nc,rs);
phi_phi = phi'*phi
phi_F = phi'*F; phi_Rs = phi'*RRs;
Ky = [1 0 0]*inv(phi_phi + R)*phi_Rs;
Kmpc = [1 0 0]*inv(phi_phi + R)*phi_F
DeltaU = inv(phi_phi + R)*phi'*(RRs - F*x)
eigenvalues = eig(A - B*Kmpc)
condition_number = cond(phi_phi + R)

%% if Np = 200
Np = 200;
disp(['Np = ', num2str(Np)]);
[F,phi,RRs] = IEmpcgains2(A,B,C,Np,Nc,rs);
phi_phi = phi'*phi
phi_F = phi'*F; phi_Rs = phi'*RRs;
Ky = [1 0 0]*inv(phi_phi + R)*phi_Rs;
Kmpc = [1 0 0]*inv(phi_phi + R)*phi_F
DeltaU = inv(phi_phi + R)*phi'*(RRs - F*x)
eigenvalues = eig(A - B*Kmpc)
condition_number = cond(phi_phi + R)

%% Note
% A large condition number means that the matrix is 
% close to being singular or non-invertible.