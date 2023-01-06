function [F, phi]= F_phi(A, B, C, Np, Nc)
F= []; phi= [];

for k= 1:Np
	phi_row= [];
	for j = 1:Nc
		if k < j
			phi_row= [phi_row 0];
		else
			phi_row = [phi_row C*A^(k-j)*B];
		end
	end
	phi = [phi; phi_row];
	F= [F; C*A^k];
end
% phi, F