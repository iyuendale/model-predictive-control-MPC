num = 0.1; den = [1 -1.4 0.48 0 0 0];
sysTF = tf(num, den)
[Am, Bm, Cm, Dm] = tf2ss(num, den)
% %%
% Am = [1.4 -0.48 0 0 0;
% 		1 0 0 0 0; 0 1 0 0 0;
% 		0 0 1 0 0; 0 0 0 1 0]
% % Am = [0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 0 -0.48 1.4]
% Bm = [1 0 0 0 0]', Cm = [0 0 0 0 0.1]
% Dm = 0;
%  augmented system
A = [Am zeros(size(Am, 1), size(Cm, 1));
	  Cm*Am eye(size(Cm, 1))]
B = [Bm; Cm*Bm], C = [zeros(size(Cm, 1), size(Am, 1)) 1]
%%  controller design
Np = 16; Nc = 4; rs = 1; rw = 0.01;
Rs = rs*ones(Np, 1); R = rw*eye(Nc);
F = []; phi = [];
for k = 1:Np
	phi_col = [];
	for i = 1:Nc
		if k < j
			phi_col = [phi_col zeros(size(C, 1), size(B, 2))];
		else
			phi_col = [phi_col C*A^(k-j)*B];
		end
	end
	phi = [phi; phi_col];
	F = [F; C*A^k];
end
F, phi, Rs, R