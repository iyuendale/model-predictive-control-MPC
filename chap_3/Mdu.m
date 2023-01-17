%% tutorial 3.7
% generates the matrix M that is applied on the difference of 
% control variable, ∆u(ki)
% M is data matrix that is applied on the first Nc samples on
% ∆u(ki)
% the block matrix Lzerot is used for constructing the control
% variable
% Nc = number of future samples for constraints to be imposed
function [M, Lzerot] = Mdu(a, N, n_in, Nc)
N_pa = sum(N);
M = zeros(n_in, N_pa);
M_du1 = zeros(n_in, N_pa);
k0 = 1;
[Al, L0] = lagd(a(k0), N(k0));
M_du1(1, 1:N(1)) = L0';
cc = N(1);

for k0 = 2:n_in
	[Al, L0] = lagd(a(k0), N(k0));
	M_du1(k0, cc+1:cc+N(k0)) = L0';
	cc = cc+N(k0);
end
Lzerot = M_du1;
M = M_du1;

for kk=2:Nc
	k0=1;
	[Al,L0]=lagd(a(k0),N(k0));
	L=Al^(kk-1)*L0;
	M_du1(1,1:N(1))=L';
	cc=N(1);
	for k0=2:n_in
		[Al,L0]=lagd(a(k0),N(k0));
		L=Al^(kk-1)*L0;
		M_du1(k0,cc+1:cc+N(k0))=L';
		cc=cc+N(k0);
	end
	M=[M;M_du1];
end