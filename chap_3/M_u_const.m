function Mu = M_u_const(a, N, Nc)
% developes matrix for inequality constraints applied on 'u'
% Nc = nomber of constraints imposed on 'u'
n_in = size(N, 1);
% M_u1 = 
start = 0; ending = 0;
M_u = zeros(Nc*n_in, sum(N));
for k = 1:n_in
	[Al, L0] = lagd(a(k), N(k));
	M_u1 = zeros(Nc, N(k));
	M_u1(1, :) = L0';
	for i = 2:Nc
		L = Al^(i-1)*L0;
		M_u1(i, :) = M_u1(i-1, :)+L';
	end
	start = 1 + ending; ending = k*N(k);
	M_u((1+(k-1)*Nc):(k*Nc), start:ending) = M_u1;
end
Mu = [M_u; -M_u];