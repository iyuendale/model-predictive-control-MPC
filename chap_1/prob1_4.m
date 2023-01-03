A = [0.6 0; 0.6 1]; B= [0.3 0.3]';C = [0 1];
x0 = [0.1 0.2]'; Nc = 4; Np = 10;
%% 1
% Y = F*x(0) + phi*deltaU
F = [];
for k = 1:Np
	F = [F; C*A^k];
end
F
phi = [];
for i = 1:Np
	phi_col = [];
	for j = 1:Nc
		if i < j
			phi_col = [phi_col zeros(size(C,1),size(B,2))];
		else
			phi_col = [phi_col C*A^(i-1)*B];
		end
	end
	phi = [phi; phi_col];
end
phi

%% 2
rs = 0; Rs = rs*ones(Np, 1); R = 3*eye(Nc, Nc);
deltaU = inv(phi'*phi + R)*phi'*(Rs - F*x0)
Y = F*x0 + phi*deltaU;
J = Y'*Y + deltaU'*R*deltaU, Y(end)

%% 3
rs = 0; Rs = rs*ones(Np, 1); R = 0*eye(Nc, Nc);
deltaU = inv(phi'*phi + R)*phi'*(Rs - F*x0)
Y = F*x0 + phi*deltaU;
J = Y'*Y + deltaU'*R*deltaU, Y(end)